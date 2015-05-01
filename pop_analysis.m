global directory slash;

%% settings
[directory,slash,user,dbldir,mapdr,servrep,mapddataf]=SetUserDir;

try
    CCNdb = connect2DB('vp_sldata');
    %     query = 'SELECT FileName FROM b_dentate';
    %     results = fetch(CCNdb,query);
    dentatefiles =fetch(CCNdb,'select r.a_file FROM recordings r WHERE r.task=''gapstop'' AND r.recloc=''dentate'''); %dentate %top_cortex
catch db_fail
    results = [];
end

%% cut file in two: make it pop_analysis_process


% CmdFileName={'S113L4A5_13500';'S114L4A5_14321';'R132L4P4_20152';'S112l4a5_12971';...
%     'S117L4A6_12741';'S118L4A5_13081';'S115L4A6_12871';...
%     'H56L5A5_21502';'S116L4A6_15450';'H53L5A5_20901';...
%     'S123L4A6_13691';'S125L4A6_13990';'H56L5A5_21321';...
%     'R97L7A1_19001'};

%% concatenate file lists if needed
% dentatefiles=[results;CmdFileName];
% dentatefiles=unique(dentatefiles);

%% prealloc
alldata=struct('fname',{},'task',{},'aligntype',{},'prevssd',{},'allmssrt_tacho',{},...
    'ssds',{},'sacdelay',{},'prefdiridx',{},...
    'pk',struct('sac',{},'vis',{},'corsac',{},'rew',{}),...
    'trialidx',struct('stoptrials',{},'goodsac',{},'noncancel',{},'cancel',{}),...
    'ndata',struct('rast',{},'alignt',{},'trialnb',{},'evttime',{}),...
    'db',struct('rec_id',{},'sort_id',{},'unit_id',{}),...
    'stats',struct('hval',{},'pval',{},'sign',{}));
%allmssrt_tacho=NaN(length(dentatefiles),1);

noclusfiles=[];

%% process files
for flbn=1:length(dentatefiles)
    dfile=dentatefiles{flbn}; %dfile=[dfile '_REX'];
    if strcmp('A',dfile(end))
        dfile=dfile(1:end-1);
    end
    
    %% subject and procdir
    if strcmp('R',dfile(1))
        subject='Rigel';
        procdir = [directory,'processed',slash,'Rigel',slash];
    elseif strcmp('S',dfile(1))
        subject='Sixx';
        procdir = [directory,'processed',slash,'Sixx',slash];
    elseif strcmp('H',dfile(1))
        subject='Hilda';
        procdir = [directory,'processed',slash,'Hilda',slash];
    end
    
    %% what file to load
    procdirlisting=dir(procdir);
    procdirfileNames={procdirlisting.name};
    loadfile=procdirfileNames(~cellfun('isempty',regexpi(procdirfileNames,dfile,'match')));
    
    % if two versions (REX and Spike2), chose REX (for now)
    if length(loadfile)>1
        loadfile=loadfile(~cellfun('isempty',regexpi(loadfile,'REX','match')));
        if isempty(loadfile) %actually it's an old file that sneaked in
            loadfile=procdirfileNames(~cellfun('isempty',regexpi(procdirfileNames,dfile,'match')));
            loadfile=loadfile(~cellfun('isempty',regexpi(loadfile,'REX','match')));
        end
    end
    try
        rectype=regexp(loadfile{:},'(?<=\d\d_)\w+(?=\.\w+)','match'); % get recording type (REX or Sp2)
    catch
        loadfile{:}
        continue
    end
    if strcmp(rectype{:},'REX')
        rectype='Rex';
    elseif strcmp(rectype{:},'Sp2')
        rectype='Spike2';
    end
    alldata(flbn,1).fname=loadfile{:}(1:end-4);
    
    %% get task and id
    try
        % if no issue with db
        query = ['SELECT r.task, r.recording_id FROM recordings r WHERE r.a_file = ''' dfile 'A'''];
        task_rec_id=fetch(CCNdb,query);
        if isempty(task_rec_id)
            %retry. Sometimes db access is deficient
            task_rec_id=fetch(CCNdb,query);
        end
        [alldata(flbn,1).task,task]=deal(task_rec_id{1});
        alldata(flbn,1).db.rec_id=task_rec_id{2};
        %check valid task
        fcodes =load([dfile,'_REX.mat'], 'allcodes');
        if ~strcmp(task,taskdetect(fcodes.allcodes))
            [alldata(flbn,1).db.rec_id,alldata(flbn,1).db.sort_id,alldata(flbn,1).db.unit_id]=deal([]);
            continue
        else
            query = ['SELECT s.sort_id FROM sorts s WHERE s.recording_id_fk = '...
                num2str(task_rec_id{2}) ' AND s.origin = ''' rectype ''' AND s.user = ''' user ''''];
            sort_id=fetch(CCNdb,query);
            if length(sort_id)==1
                usesort=1;
                alldata(flbn,1).db.sort_id=sort_id{1};
            else
                %% need to choose which sort!
                disp('too many sorts to choose from')
                %             usesort=1;
                break
            end
            % get the unit ID (cluster here means spike cluster - from spike sorting)
            query = ['SELECT c.cluster_id FROM clusters c WHERE c.sort_id_fk = '...
                num2str(sort_id{usesort}) ' AND c.phenotype != ''junk'''];
            unit_id=fetch(CCNdb,query);
            if isempty(unit_id)
                % create one, e basta
                [~, alldata(flbn,1).db.unit_id] = addCluster(sort_id{usesort}, 1,CCNdb);
%                 noclusfiles=[noclusfiles; dfile];
                continue
            elseif length(unit_id)==1
                useunit=1;
                alldata(flbn,1).db.unit_id=unit_id{1};
            elseif length(unit_id)> 1%need to go 1 by 1
                disp('need to process units 1 by 1')
                break
            end
        end
        clear fcodes;
    catch db_fail
        task='gapstop';
        r_id=[];
    end
    %% sort ou
    if strcmp(task,'st_saccades')
        %test peak shift on prefered dir vs anti-dir
    elseif strcmp(task,'gapstop')
        
        %% align rasters
        % common presets
        [~, trialdirs] = data_info(loadfile{:}, 1, 1); %reload file: yes (shouldn't happen, though), skip unprocessed files: yes
        includebad=0;
        spikechannel=1; %select appropriate cluster
        keepdir='compall'; %alldir %which sac directions
        togrey=[];
        singlerastplot=0;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% prealloc
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fails=[];
        alignments={'failed_fast','correct_slow','ssd','corrsacfailed','rewcorrect_rewslow'};
        for algn=1:5
            %% set parameter values
            if algn==1 % sac vs stop: 3 alignements 'sac' (correct sac) / 'stop_cancel' (to SS + SSRT) / 'stop_non_cancel' (incorrect sac)
                firstalign=6;
                secondalign=8;
                plottype = 0;
                plotstart=1100;
                plotstop=1000;
                option=NaN;
            elseif algn==2 % tgt vs stop
                firstalign=7;
                secondalign=8;
                plottype = 0; % 3 for splitting data in three groups: short SSD, med SSD and long SSD
                plotstart=200;
                plotstop=600;
                if plottype == 3;
                    option='truealign';
                else
                    option=NaN;
                end
            elseif algn==3 % ssd
                firstalign=7; % as if align to target
                secondalign=507;
                plottype = 0;
                plotstart=800;
                plotstop=800;
                
                %% ssd alignement specific:
                % if aligning to ssd, need to align NSS trials according to latency
                
                % get psychometric values
                [mssrt,inhibfun,ccssd,nccssd,ssdvalues,tachomc,tachowidth,...
                    sacdelay,rewtimes,prevssd,trialidx]=findssrt(loadfile{:}, 0);
                if ~isnan(mssrt)
                    mssrt=max([mssrt (mean(tachomc)+tachowidth/2)]); %replace by: if mssrt < tachomc+tachowidth/2, mssrt=tachomc+tachowidth/2, end; ?
                    if mssrt> 130 && mean(tachomc)>50
                        mssrt=mean(tachomc)+tachowidth;
                    end
                else
                    alldata(flbn,1).allmssrt_tacho=NaN;
                    continue
                end
                alldata(flbn,1).allmssrt_tacho={mssrt mean(tachomc) tachowidth};
                alldata(flbn,1).prevssd={prevssd};
                alldata(flbn,1).sacdelay={sacdelay};
                alldata(flbn,1).trialidx=trialidx;
                
                if size(trialidx.goodsac,1)~=size(alldata(flbn, 1).ndata(1, 1).rast ,1)
                    if size(trialidx.goodsac,1)>size(alldata(flbn, 1).ndata(1, 1).rast,1)
                        alldata(flbn,1).sacdelay={alldata(flbn,1).sacdelay{:}(ismember(trialidx.goodsac,alldata(flbn,1).ndata(1,1).trialnb))};
                    elseif size(trialidx.goodsac,1)<size(alldata(flbn, 1).ndata(1, 1).rast,1) %that may be an issue later
                        %let's see if it happens, first
                        size(trialidx.goodsac,1)
                    end
                end
                
                %allssds=unique([ccssd;nccssd]);
                
                % canceled trials
                ccssdval=unique(ccssd);
                ctmatchlatidx=zeros(length(sacdelay),length(ccssdval));
                for ssdval=1:length(ccssdval)
                    % Keeping NSS trials with sac latencies long enough
                    % that they would have occured after a stop-signal
                    ctmatchlatidx(:,ssdval)=sacdelay>ccssdval(ssdval)+round(mssrt);
                end
                nullidx=sum(ctmatchlatidx,2)==0;
                ctmatchlatidx(nullidx,1)=1;
                % getting ssds for each NNS trial, taking the highest ssd.
                ctmatchlatidx=ccssdval(sum(ctmatchlatidx,2));
                ctmatchlatidx(nullidx,1)=0;
                
                % non-canceled trials
                nccssdval=sort(unique(nccssd));
                nctallmatchlatidx=zeros(length(sacdelay),length(nccssdval));
                for ssdval=1:length(nccssdval) % keeping NSS trials in which
                    % a saccade would have been
                    % initiated even if a stop
                    % signal had occurred, but
                    % with saccade latencies
                    % greater than the
                    % stop-signal delay plus a
                    % visual-response latency.
                    % We take tachomc-tachowidth/2
                    % rather than the arbitrary
                    % 50ms from Hanes et al 98
                    nctallmatchlatidx(:,ssdval)=sacdelay>nccssdval(ssdval)+(mean(tachomc)-tachowidth/2) & sacdelay<nccssdval(ssdval)+round(mssrt);
                end
                % getting ssds for each NNS trial, taking the lowest ssd.
                nctmatchlatidx=zeros(size(nctallmatchlatidx,1),1);
                for midx=1:size(nctallmatchlatidx,1)
                    if ~isempty(find(nctallmatchlatidx(midx,:),1))
                        nctmatchlatidx(midx)=nccssdval(find(nctallmatchlatidx(midx,:),1));
                    end
                end
                option=[ctmatchlatidx nctmatchlatidx];
                
            elseif algn==4 % align to corrective saccades in failed cancellation trial
                firstalign=9; % corrective saccades for SST
                secondalign=8; % will find corrective saccades in failed cancellation trial
                plottype = 0;
                plotstart=1000;
                plotstop=500;
                option=NaN;
            elseif algn==5 % align to reward time for NSS and CS trials
                firstalign=4;
                % will find reward time in cancelled saccade
                % for failed cancellation, two possibility:
                % align to error code (easy, but not time consistent with other trial types)
                % align to estimated reward time (time between saccade and
                % reward, calculated from NSS trials)
                secondalign=8;
                plottype = 0;
                plotstart=1000;
                plotstop=200;
                option=NaN;
            end
            
            alldata(flbn,algn).aligntype=alignments{algn};
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% align rasters to conditions
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %use GUI-independent prealign
            getaligndata={}; %re-init structure
            try
                getaligndata = prealign(loadfile{:}(1:end-4), trialdirs, task, firstalign,...
                    secondalign,  includebad, spikechannel, keepdir,...
                    togrey, singlerastplot, option); % align data, don't plot rasters
            catch prealign_fail
                fails={fails; [loadfile{:}(1:end-4), prealign_fail.message]}; %prealign_fail
                continue
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% z-score pre-ssd and pre-sac?
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% get peak firing rate for future normalization. Find prefered dir.
            % also keep ssds for cs and ncs trials
            if find(strcmp({getaligndata.alignlabel},'sac'))
                numrastrow=arrayfun(@(x) size(x.rasters,1), getaligndata, 'UniformOutput', false);
                colrast=nan(sum([numrastrow{:}]),601);
                colrastidx=[0 numrastrow{:}];
                prefdir=nan(2,3);
                for alignd=1:2:size(getaligndata,2) %align on trials where saccades occured
                    if colrastidx(alignd+1)==0
                        continue
                    end
                    sacalgrasters=getaligndata(1,alignd).rasters;
                    alignmtt=getaligndata(1,alignd).alignidx;
                    start=alignmtt-300; stop=alignmtt+300; % -300 to 300 time window around sac (at 0).
                    colrast(colrastidx(alignd)+sum(colrastidx(1:alignd-1))+1:colrastidx(alignd+1)+sum(colrastidx(1:alignd)),:)=...
                        sacalgrasters(:,start:stop);
                    % prefered dir
                    [unikdir,~,unikidx]=unique(getaligndata(1,alignd).dir);
                    convpkdir=nan(length(unikdir),1);
                    for convdirs=1:length(unikdir)
                        convpkdir(convdirs,:)=max(conv_raster(sacalgrasters(unikidx==convdirs,start:stop),20));
                    end
                    if ~isempty(convpkdir)
                        if length(unikdir(round(convpkdir)==max(round(convpkdir))))>2
                            prefdir(1:2,alignd)=unikdir((convpkdir)==max((convpkdir)));
                        else
                            prefdir(1:2,alignd)=unikdir(round(convpkdir)==max(round(convpkdir))); %if no prefered dir, will be decided just below
                        end
                    else
                        prefdir(alignd)=NaN;
                    end
                end
                %get most found pref dire
                [~,~,prefdir]=mode(prefdir(~isnan(prefdir)));
                if length(prefdir{:})>1 %then keep dir with most trials
                    prefdirnbtrial=([sum([getaligndata.dir]==prefdir{:}(1)) sum([getaligndata.dir]==prefdir{:}(2))]);
                    prefdir={prefdir{:}(find(prefdirnbtrial==max(prefdirnbtrial),1))};
                end
                alldata(flbn,algn).prefdiridx=arrayfun(@(x) x.dir==prefdir{:},getaligndata,'UniformOutput',false);
                convrasters=conv_raster(colrast,10,1,size(colrast,2));
                alldata(flbn,algn).pk.sac=max(convrasters);
                pk_or_tro_time=find(convrasters==max(convrasters) | convrasters==min(convrasters),1);
                %% make some stats on sac alignment
                try
                    sacalgrasters=getaligndata(1,1).rasters;
                    alignmtt=getaligndata(1,1).alignidx;
                    start=alignmtt-600+pk_or_tro_time; stop=alignmtt+pk_or_tro_time;
                    % test statcond pre/post peak_or_trough, and direction of change.
                    [alldata(flbn,algn).stats.hval, alldata(flbn,algn).stats.pval,...
                        alldata(flbn,algn).stats.sign]=rastplotstat(sacalgrasters,10,...
                        [alignmtt-600+pk_or_tro_time alignmtt-(300-pk_or_tro_time)+30],...
                        [alignmtt-(300-pk_or_tro_time)+30 alignmtt+pk_or_tro_time],0);
                catch
                    [alldata(flbn,algn).stats.hval, alldata(flbn,algn).stats.pval,...
                        alldata(flbn,algn).stats.sign]=deal(NaN);
                end
                
            elseif find(strcmp({getaligndata.alignlabel},'corsac'))
                if size(getaligndata,2)>2 && size(getaligndata(3).rasters,1)>1
                    numrastrow=arrayfun(@(x) size(x.rasters,1), getaligndata, 'UniformOutput', false);
                    colrast=nan(sum([numrastrow{:}]),601);
                    colrastidx=[0 numrastrow{:}];
                    prefdir=nan(2,3);
                    for alignd=1:2:3
                        if colrastidx(alignd+1)==0
                            continue
                        end
                        sacalgrasters=getaligndata(1,alignd).rasters;
                        alignmtt=getaligndata(1,alignd).alignidx;
                        start=alignmtt-300; stop=alignmtt+300; % -300 to 300 time window around sac (at 0).
                        colrast(colrastidx(alignd)+sum(colrastidx(1:alignd-1))+1:colrastidx(alignd+1)+sum(colrastidx(1:alignd)),:)=...
                            sacalgrasters(:,start:stop);
                        % prefered dir
                        [unikdir,~,unikidx]=unique(getaligndata(1,alignd).dir);
                        convpkdir=nan(length(unikdir),1);
                        for convdirs=1:length(unikdir)
                            convpkdir(convdirs,:)=max(conv_raster(sacalgrasters(unikidx==convdirs,start:stop),20));
                        end
                        if ~isempty(convpkdir)
                            if length(unikdir(round(convpkdir)==max(round(convpkdir))))>2
                                prefdir(1:2,alignd)=unikdir((convpkdir)==max((convpkdir)));
                            else
                                prefdir(1:2,alignd)=unikdir(round(convpkdir)==max(round(convpkdir))); %if no prefered dir, will be decided just below
                            end
                        else
                            prefdir(alignd)=NaN;
                        end
                    end
                    %get most found pref dir
                    [~,~,prefdir]=mode(prefdir(~isnan(prefdir)));
                    if length(prefdir{:})>1 %then keep dir with most trials
                        prefdirnbtrial=([sum([getaligndata.dir]==prefdir{:}(1)) sum([getaligndata.dir]==prefdir{:}(2))]);
                        prefdir={prefdir{:}(find(prefdirnbtrial==max(prefdirnbtrial),1))};
                    end
                    alldata(flbn,algn).prefdiridx=arrayfun(@(x) x.dir==prefdir{:},getaligndata,'UniformOutput',false);
                    convrasters=conv_raster(colrast,10,11,size(colrast,2)-10);
                    alldata(flbn,algn).pk.corsac=max(convrasters);
                end
            elseif find(strcmp({getaligndata.alignlabel},'tgt'))
                numrastrow=arrayfun(@(x) size(x.rasters,1), getaligndata, 'UniformOutput', false);
                colrast=nan(sum([numrastrow{:}]),251);
                colrastidx=[0 numrastrow{:}];
                prefdir=nan(2,size(getaligndata,2));
                for alignd=1:size(getaligndata,2)
                    if colrastidx(alignd+1)==0
                        continue
                    end
                    sacalgrasters=getaligndata(1,alignd).rasters;
                    alignmtt=getaligndata(1,alignd).alignidx;
                    start=alignmtt; stop=alignmtt+250; % -300 to 300 time window around sac (at 0).
                    colrast(colrastidx(alignd)+sum(colrastidx(1:alignd-1))+1:colrastidx(alignd+1)+sum(colrastidx(1:alignd)),:)=...
                        sacalgrasters(:,start:stop);
                    % prefered dir
                    [unikdir,~,unikidx]=unique(getaligndata(1,alignd).dir);
                    convpkdir=nan(length(unikdir),1);
                    for convdirs=1:length(unikdir)
                        convpkdir(convdirs,:)=max(conv_raster(sacalgrasters(unikidx==convdirs,start:stop),20));
                    end
                    if ~isempty(convpkdir)
                        if length(unikdir(round(convpkdir)==max(round(convpkdir))))>2
                            prefdir(1:2,alignd)=unikdir((convpkdir)==max((convpkdir)));
                        else
                            prefdir(1:2,alignd)=unikdir(round(convpkdir)==max(round(convpkdir))); %if no prefered dir, will be decided just below
                        end
                    else
                        prefdir(alignd)=NaN;
                    end
                end
                %get most found pref dir
                [~,~,prefdir]=mode(prefdir(~isnan(prefdir)));
                if length(prefdir{:})>1 %then keep dir with most trials
                    prefdirnbtrial=([sum([getaligndata.dir]==prefdir{:}(1)) sum([getaligndata.dir]==prefdir{:}(2))]);
                    prefdir={prefdir{:}(find(prefdirnbtrial==max(prefdirnbtrial),1))};
                end
                alldata(flbn,algn).prefdiridx=arrayfun(@(x) x.dir==prefdir{:},getaligndata,'UniformOutput',false);
                convrasters=conv_raster(colrast,10,11,size(colrast,2)-10);
                alldata(flbn,algn).pk.vis=max(convrasters);
                % keep ssds
                stoptrialsdat=cellfun(@(x) strfind(x,'stop'),{getaligndata(1,:).alignlabel},'UniformOutput',false);
                stoptrialsdat=getaligndata(1,~cellfun('isempty',stoptrialsdat));
                if isfield(stoptrialsdat,'ssd') && ~isempty(cat(1,stoptrialsdat.ssd))
                    alldata(flbn,1).ssds={stoptrialsdat.ssd};
                else
                    alldata(flbn,1).ssds=[];
                end
            elseif find(strcmp({getaligndata.alignlabel},'rew'))
                numrastrow=arrayfun(@(x) size(x.rasters,1), getaligndata, 'UniformOutput', false);
                colrast=nan(sum([numrastrow{:}]),501);
                colrastidx=[0 numrastrow{:}];
                prefdir=nan(2,2);
                for alignd=1:2
                    if colrastidx(alignd+1)==0
                        continue
                    end
                    sacalgrasters=getaligndata(1,alignd).rasters;
                    alignmtt=getaligndata(1,alignd).alignidx;
                    start=alignmtt-300; stop=alignmtt+200; % -300 to 300 time window around sac (at 0).
                    colrast(colrastidx(alignd)+sum(colrastidx(1:alignd-1))+1:colrastidx(alignd+1)+sum(colrastidx(1:alignd)),:)=...
                        sacalgrasters(:,start:stop);
                    % prefered dir
                    [unikdir,~,unikidx]=unique(getaligndata(1,alignd).dir);
                    convpkdir=nan(length(unikdir),1);
                    for convdirs=1:length(unikdir)
                        convpkdir(convdirs,:)=max(conv_raster(sacalgrasters(unikidx==convdirs,start:stop),20));
                    end
                    if ~isempty(convpkdir)
                        if length(unikdir(round(convpkdir)==max(round(convpkdir))))>2
                            prefdir(1:2,alignd)=unikdir((convpkdir)==max((convpkdir)));
                        else
                            prefdir(1:2,alignd)=unikdir(round(convpkdir)==max(round(convpkdir))); %if no prefered dir, will be decided just below
                        end
                    else
                        prefdir(alignd)=NaN;
                    end
                end
                %get most found pref dir
                [~,~,prefdir]=mode(prefdir(~isnan(prefdir)));
                if length(prefdir{:})>1 %then keep dir with most trials
                    prefdirnbtrial=([sum([getaligndata.dir]==prefdir{:}(1)) sum([getaligndata.dir]==prefdir{:}(2))]);
                    prefdir={prefdir{:}(find(prefdirnbtrial==max(prefdirnbtrial),1))};
                end
                alldata(flbn,algn).prefdiridx=arrayfun(@(x) x.dir==prefdir{:},getaligndata,'UniformOutput',false);
                convrasters=conv_raster(colrast,10,11,size(colrast,2)-10);
                alldata(flbn,algn).pk.rew=max(convrasters);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%
            %% store rasters, alignment index, trial index, "grey" event times
            [alldata(flbn,algn).ndata(1:size({getaligndata.rasters},2)).rast]=deal(getaligndata.rasters);
            [alldata(flbn,algn).ndata(1:size({getaligndata.rasters},2)).alignt]=deal(getaligndata.alignidx);
            [alldata(flbn,algn).ndata(1:size({getaligndata.rasters},2)).trialnb]=deal(getaligndata.trials);
            [alldata(flbn,algn).ndata(1:size({getaligndata.rasters},2)).evttime]=deal(getaligndata.allgreyareas);
            
            % [t df pvals] = statcond({convrasters closeconvrasters}, 'method', 'perm', 'naccu', 20000,'verbose','off');
            
            %% store data for population plotting
            
            
        end
        
        %% classification: ramp/burst vs pause rebound
        %% Normalization by peak firing
    elseif strcmp(task,'tokens')
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% analyze gapstop data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gsdlist=cellfun(@(x) strcmp(x,'gapstop'),{alldata(:,1).task}) & ~cellfun('isempty',{alldata(:,1).ndata});
recluster=0;

%% reshape data
gsdata.allgsalignmnt=reshape({alldata(gsdlist,:).aligntype},size(alldata(gsdlist,:)));
gsdata.allgsprevssd=reshape({alldata(gsdlist,:).prevssd},size(alldata(gsdlist,:)));gsdata.allgsprevssd=gsdata.allgsprevssd(:,1);
gsdata.allgsssds=reshape({alldata(gsdlist,:).ssds},size(alldata(gsdlist,:)));gsdata.allgsssds=gsdata.allgsssds(:,1);
gsdata.allgssacdelay=reshape({alldata(gsdlist,:).sacdelay},size(alldata(gsdlist,:)));gsdata.allgssacdelay=gsdata.allgssacdelay(:,1);
gsdata.allgsprefdir=reshape({alldata(gsdlist,:).prefdiridx},size(alldata(gsdlist,:)));
gsdata.allgsndata=reshape({alldata(gsdlist,:).ndata},size(alldata(gsdlist,:)));
gsdata.allgsmssrt_tacho=reshape({alldata(gsdlist,:).allmssrt_tacho},size(alldata(gsdlist,:))); gsdata.allgsmssrt_tacho=gsdata.allgsmssrt_tacho(:,1);
gsdata.alldb=reshape({alldata(gsdlist,:).db},size(alldata(gsdlist,:))); gsdata.alldb=gsdata.alldb(:,1);

% allgstasks=reshape({alldata(gsdlist,:).task},size(alldata(gsdlist,:))); allgstasks=allgstasks(:,1);
% allgspk=reshape({alldata(gsdlist,:).pk},size(alldata(gsdlist,:)));
% allgs_rec_id=reshape({alldata(gsdlist,:).db_rec_id},size(alldata(gsdlist,:)));
% allgsstats=reshape({alldata(gsdlist,:).stats},size(alldata(gsdlist,:)));
% allgstrialidx=reshape({alldata(gsdlist,:).trialidx},size(alldata(gsdlist,:)));
% allgsfname=reshape({alldata(gsdlist,:).fname},size(alldata(gsdlist,:)));

% load('E:\BoxSync\Box Sync\Home Folder vp35\Sync\SommerLab\projects\countermanding\popclusters\countermanding_cDn_gsdata.mat')

pop_a_countermanding(gsdata,recluster,CCNdb);

