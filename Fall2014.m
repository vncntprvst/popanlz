
%% settings
[directory,slash,user,dbldir,mapdr,servrep,mapddataf]=SetUserDir;

CCNdb = connect2DB('vp_sldata');

query = 'SELECT FileName FROM b_dentate';
results = fetch(CCNdb,query);

CmdFileName={'S113L4A5_13500';'S114L4A5_14321';'R132L4P4_20152';'S112l4a5_12971';...
    'S117L4A6_12741';'S118L4A5_13081';'S115L4A6_12871';...
    'H56L5A5_21502';'S116L4A6_15450';'H53L5A5_20901';...
    'S123L4A6_13691';'S125L4A6_13990';'H56L5A5_21321';...
    'R97L7A1_19001'};

%% concatenate file lists
dentatefiles=[results;CmdFileName];
dentatefiles=unique(dentatefiles);

%% prealloc
alldata=struct('task',{},'aligntype',{},'allmssrt',{},...
    'pk',struct('sac',{},'vis',{}),'ndata',{});
%allmssrt=NaN(length(dentatefiles),1);


%% process files
for flbn=1:length(dentatefiles)
    dfile=dentatefiles{flbn}; %dfile=[dfile '_REX'];
    
    %% get task and id
    query = ['SELECT r.task, r.recording_id FROM recordings r WHERE r.a_file = ''' dfile 'A'''];
    results=fetch(CCNdb,query);
    task=results{1};
    r_id=results{2};
    
    %% sort ou
    if strcmp(task,'st_saccades')
        %test peak shift on prefered dir vs anti-dir
    elseif strcmp(task,'gapstop')
        alldata(flbn,1).task=task;
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
        
        % if two versions (REX and Spike2), chose REX.
        if length(loadfile)>1
            loadfile=loadfile(~cellfun('isempty',regexpi(loadfile,'REX','match')));
            if isempty(loadfile) %actually it's an old file that sneaked in
                loadfile=procdirfileNames(~cellfun('isempty',regexpi(procdirfileNames,dfile,'match')));
                loadfile=loadfile(~cellfun('isempty',regexpi(loadfile,'REX','match')));
            end
        end
        
        %% first, get psychometric values
        [mssrt,inhibfun,ccssd,nccssd,ssdvalues,tachomc,tachowidth,sacdelay,rewtimes]=findssrt(loadfile{:}, 0);
        mssrt=max([mssrt (mean(tachomc)+tachowidth/2)]); %replace by: if mssrt < tachomc+tachowidth/2, mssrt=tachomc+tachowidth/2, end; ?
        %         if mssrt>(mean(tachomc)+3*tachowidth)
        %             mssrt=mean(tachomc)+tachowidth;
        %         end
        alldata(flbn,1).allmssrt=mssrt;
        
        %% align rasters
        % common presets
        [~, trialdirs] = data_info(loadfile{:}, 1, 1); %reload file: yes (shouldn't happen, though), skip unprocessed files: yes
        includebad=0;
        spikechannel=1; %select appropriate cluster
        keepdir='compall'; %alldir %which sac directions
        togrey=[];
        singlerastplot=0;
        
        %%prealloc
                
        for alignment=1:3
            %% set parameter values
            if alignment==1 % sac vs stop
                firstalign=6;
                secondalign=8;
                aligntype='failed_fast';
                plottype = 0;
                plotstart=1000;
                plotstop=1000;
                option=NaN;
            elseif alignment==2 % tgt vs stop
                firstalign=7;
                secondalign=8;
                aligntype='correct_slow';
                plottype = 0; % 3 for splitting data in three groups: short SSD, med SSD and long SSD
                plotstart=200;
                plotstop=600;
                if plottype == 3;
                    option='truealign';
                else
                    option=NaN;
                end
            elseif alignment==3 % ssd
                firstalign=7; % as if align to target
                secondalign=507;
                aligntype='ssd';
                plottype = 0;
                plotstart=800;
                plotstop=600;
                
                %% ssd alignement specific:
                % if aligning to ssd, got to align NSS trials according to latency
                %allssds=unique([ccssd;nccssd]);
                
                % canceled trials
                ccssdval=unique(ccssd);
                ctmatchlatidx=zeros(length(sacdelay),length(ccssdval));
                for ssdval=1:length(ccssdval)
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
                for ssdval=1:length(nccssdval)
                    nctallmatchlatidx(:,ssdval)=sacdelay>nccssdval(ssdval)+50 & sacdelay<nccssdval(ssdval)+round(mssrt);
                end
                % getting ssds for each NNS trial, taking the lowest ssd.
                nctmatchlatidx=zeros(size(nctallmatchlatidx,1),1);
                for midx=1:size(nctallmatchlatidx,1)
                    if ~isempty(find(nctallmatchlatidx(midx,:),1))
                        nctmatchlatidx(midx)=nccssdval(find(nctallmatchlatidx(midx,:),1));
                    end
                end
                option=[ctmatchlatidx nctmatchlatidx];
            end
            
            alldata(flbn,alignment).aligntype=aligntype;
            
            %% use GUI-independent prealign
            try
            getaligndata = prealign(loadfile{:}(1:end-4), trialdirs, task, firstalign,...
                secondalign,  includebad, spikechannel, keepdir,...
                togrey, singlerastplot, option); % align data, don't plot rasters
            catch
                continue
            end
            %% z-score pre-ssd and pre-sac? 
            
            %% get peak firing rate for future nomrmalization
            if find(strcmp({getaligndata.alignlabel},'sac'))
                sacalgrasters=getaligndata(1,strcmp({getaligndata.alignlabel},'sac')).rasters;
                alignmtt=getaligndata(1,strcmp({getaligndata.alignlabel},'sac')).alignidx;
                start=alignmtt-300; stop=alignmtt+300; % -300 to 300 time window around sac (at 0).
                convrasters=conv_raster(sacalgrasters,10,start,stop);
                alldata(flbn,alignment).pk.sac=nanmean(convrasters);
            elseif find(strcmp({getaligndata.alignlabel},'tgt'))
                sacalgrasters=getaligndata(1,strcmp({getaligndata.alignlabel},'tgt')).rasters;
                alignmtt=getaligndata(1,strcmp({getaligndata.alignlabel},'tgt')).alignidx;
                start=alignmtt; stop=alignmtt+250; % 0 to 250 time window around cue (at 0).
                convrasters=conv_raster(sacalgrasters,10,start,stop);
                alldata(flbn,alignment).pk.vis=nanmean(convrasters);
            end
            
            alldata(flbn,alignment).ndata=getaligndata;
            
            
% [t df pvals] = statcond({convrasters closeconvrasters}, 'method', 'perm', 'naccu', 20000,'verbose','off');
 
            
             
            %% store data for population plotting 
            
            
        end
        
        %% classification: ramp/burst vs pause rebound
        %% Normalization by peak firing
    elseif strcmp(task,'tokens')
    end
    
end

alltasks=reshape({alldata.task},size(alldata));
allalignmnt=reshape({alldata.aligntype},size(alldata));
allmssrt=reshape({alldata.allmssrt},size(alldata));
allpk=reshape({alldata.pk},size(alldata));
allndata=reshape({alldata.ndata},size(alldata));





% outputs = struct('mssrt',{},...
%     'ssdvalues',{});

% datainsert(conn,'recordings',col_names, this_data);
%         commit(conn);