function pop_a_countermanding_tgtalgn_specific(data,proc_option,conn)
global directory slash;
if isempty(directory)
    [directory,slash]=SetUserDir;
end

% processing options
if proc_option.singlessd
    proc_option.prefdironly=0; %otherwise we don't keep anything
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% normalize data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1/ convolve rasters with 200ms before saccade, 200 after saccade, 20ms kernel
%time window. Add kernel * 6 ms (see fullgauss_filtconv), e.g. 60 ms at both
% ends, which will be cut.
% data.allndata has 3 column for 3 aligntype. Each cell has 3 or 4 for different conditions
sigma=10;
baslineLength=500;
[sacresps,sacrespsTrials]=cellfun(@(x) conv_raster(x(1,1).rast,sigma,x(1,1).alignt-(200+sigma*3),x(1,1).alignt+(199+sigma*3)), data.allndata(:,1), 'UniformOutput',false); %400ms period
[bslresps,bslrespsTrials]=cellfun(@(x) conv_raster(x(1,1).rast,sigma,x(1,1).alignt-(baslineLength+sigma*3),x(1,1).alignt+(sigma*3-1)), data.allndata(:,2), 'UniformOutput',false); %500ms period
fullresps=cellfun(@(x) conv_raster(x(1,1).rast,sigma,1,size(x(1,1).rast,2)), data.allndata(:,2), 'UniformOutput',false); %full response

%% remove bad apples
badapl=cellfun(@(x) size(x,2)==1, sacresps);
sacresps=sacresps(~badapl,:);
sacresps=cat(1,sacresps{:});
bslresps=bslresps(~badapl,:);
bslresps=cat(1,bslresps{:});
fullresps=fullresps(~badapl,:);
% fullresps=cat(1,fullresps{:});

% clusterIdx=clusterIdx(~badapl,:);
% figure; plot(mean(sacresps(clusterIdx==5,:)))

data.allalignmnt=data.allalignmnt(~badapl,:);
data.allmssrt_tacho=data.allmssrt_tacho(~badapl,1);
%allpk=allpk(~badapl,:); %not needed
data.allndata=data.allndata(~badapl,:);
%all_rec_id=all_rec_id(~badapl,1); %not needed
%allstats=allstats(~badapl,1); %not needed
data.allprevssd=data.allprevssd(~badapl,:);
data.allssds=data.allssds(~badapl,:);
data.allsacdelay=data.allsacdelay(~badapl,:);
data.allprefdir=data.allprefdir(~badapl,:);
%alltrialidx=alltrialidx(~badapl,:); %not needed
data.alldb=data.alldb(~badapl,:);

% 2/ standardize response
% z-score normalization by baseline - based on pre-target activity
bslresp_mean=nanmean(bslresps');
bslresp_sd=nanstd(bslresps');
% bnorm_sacresps is used for clustering purposes only
bnorm_sacresps=(sacresps-repmat(bslresp_mean',1,size(sacresps,2)))./repmat(bslresp_sd',1,size(sacresps,2));

% z-score normalization over response period (alternative method, if typically low
% baseline firing rate). Also forces clustering to operate on shapes rather than
% amplitude, by squashing response range
sacresp_mean=nanmean(sacresps');
sacresp_sd=nanstd(sacresps');
rnorm_sacresps=(sacresps-repmat(sacresp_mean',1,size(sacresps,2)))./repmat(sacresp_sd',1,size(sacresps,2));

% full response norm
fr_mean=cellfun(@(x) nanmean(x),fullresps);
fr_sd=cellfun(@(x) nanstd(x),fullresps);

%plot standardized population
% figure;
% imagesc(rnorm_sacresps);
% colormap gray
% colorbar;

%% look at "best" cells
% midrange=size(bnorm_sacresps,2)/2;
% % Make seed that represent midrange drop / midrange burst / ramp to end / ramp down
% midrangedropseeds=cellfun(@(x) mean(x(1,midrange-150:midrange-50))-mean(x(1,midrange+50:midrange+150)), mat2cell(bnorm_sacresps,ones(size(bnorm_sacresps,1),1)));
% % ramp to end
% outerrangerampseeds=cellfun(@(x) mean(x(1,length(x)-150:length(x)-1))-mean(x(1,1:150)), mat2cell(bnorm_sacresps,ones(size(bnorm_sacresps,1),1)));
% % for ramps all the way down, keep only non-bursting / falling response (~monotonic)
% leastdiff_bnorm_sacresps=bnorm_sacresps(max(abs(diff(bnorm_sacresps)),[],2)<5,:);
% outerrangerampdownseeds=cellfun(@(x) mean(x(1,length(x)-150:length(x)-1))-mean(x(1,1:150)), mat2cell(leastdiff_bnorm_sacresps,ones(size(leastdiff_bnorm_sacresps,1),1)));
% % diff sort works for peaks as well, by opposition, and could be used
% % to separate sharp bursts from smoth bursts (and template 2 from 3 apparently):
% % [~,pkseeds_vals_idx]=sort(max(abs(diff(bnorm_sacresps)),[],2),'descend');
% midrangepeakseeds=cellfun(@(x) (mean(x(1,midrange+50:midrange+100))-mean(x(1,midrange-150:midrange-50)))+...
%     (mean(x(1,midrange+50:midrange+100))-mean(x(1,midrange+100:midrange+200))), mat2cell(bnorm_sacresps,ones(size(bnorm_sacresps,1),1)));
%
% % keep 10 highest seed values
% [~,mrdropseeds_vals_idx]=sort(midrangedropseeds);
% [~,mrpkseeds_vals_idx]=sort(midrangepeakseeds);
% [~,orruseeds_vals_idx]=sort(outerrangerampseeds);
% [~,orrdseeds_vals_idx]=sort(outerrangerampdownseeds);
% top_drop=mrdropseeds_vals_idx(end-10:end);
% top_burst=mrpkseeds_vals_idx(end-10:end);
% top_rampatw=orruseeds_vals_idx(end-10:end);
% top_rampdown=orrdseeds_vals_idx(1:11);
%
% figure;
% for topfig=1:size(top_drop,1)
%     try
%     align=data.allndata{top_drop(topfig), 3}(4).alignt;
%     rasters=((data.allndata{top_drop(topfig), 3}(4).rast(:,align-800:align+800)));
%     subplot(2,1,2)
%     hold on
%     plot(conv_raster(rasters))
%     subplot(2,1,1)
%     [indy, indx] = ind2sub(size(rasters),find(rasters));
%     plot([indx';indx'],[indy';indy'+1],'LineStyle','-'); % plot rasters
%     catch
%         continue
%     end
% end


%%  get unit cluster info and profiles

unit_ids=cellfun(@(x) x.unit_id,data.alldb);


[sorted_unit_ids,sunitid_idx]=sort(unit_ids);
query = ['SELECT c.profile, c.profile_type FROM clusters c WHERE cluster_id IN (' sprintf('%.0f,' ,sorted_unit_ids(1:end-1)') num2str(sorted_unit_ids(end)) ')'];
profiles = fetch(conn,query);
sunitid_revidx(sunitid_idx)=1:length(unit_ids);
clusidx=[profiles{sunitid_revidx,2}];
clustypes={profiles{sunitid_revidx,1}};


if proc_option.ssdpkalign==1
    % Instead, classify according to when peak occurs with respect to stop signal
    % data span
    ssd_startstop=[800 700];
    
    % kernel
    conv_sigma=50;
    half_sixsig=conv_sigma*3; %half kernel window
    
    for gsd=1:size(data.allndata,1)
        if clusidx(gsd)~=-1
            try
                rasters=data.allndata{gsd,3}(3).rast;
                alignmtt=data.allndata{gsd,3}(3).alignt;
                start=alignmtt-ssd_startstop(1)-half_sixsig; stop=alignmtt+ssd_startstop(2)+half_sixsig;
                normsdf=conv_raster(rasters,conv_sigma,start,stop,bslresps(gsd,:)); %normalize by baseline
                % find peak and drough
                peak_drough=find(abs(round(diff(fullgauss_filtconv(normsdf,5,0,'same')),4).*1000)<=3);
                % segment sequence in pre, peri , post, late
                pre_ssd=normsdf(peak_drough(peak_drough>ssd_startstop(1)-600 & peak_drough<ssd_startstop(1)-200)); if isempty(pre_ssd); pre_ssd=0; end
                peri_ssd=normsdf(peak_drough(peak_drough>ssd_startstop(1)-200 & peak_drough<ssd_startstop(1)+100)); if isempty(peri_ssd); peri_ssd=0; end
                post_ssd=normsdf(peak_drough(peak_drough>ssd_startstop(1)+100 & peak_drough<ssd_startstop(1)+400)); if isempty(post_ssd); post_ssd=0; end
                late_trial=normsdf(peak_drough(peak_drough>ssd_startstop(1)+400 & peak_drough<max(normsdf))); if isempty(late_trial); late_trial=0; end
                
                if max(peri_ssd) > max([pre_ssd post_ssd]) % Clus 1
                    clusidx(gsd)=101;
                elseif max(post_ssd) > max([pre_ssd peri_ssd late_trial]) % Clus 2
                    clusidx(gsd)=102;
                elseif max(late_trial) > max([pre_ssd peri_ssd post_ssd]) % Clus 3
                    clusidx(gsd)=103;
                elseif max(pre_ssd) > max([peri_ssd post_ssd late_trial]) % Clus 4
                    clusidx(gsd)=104;
                end
            catch
                continue
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Separate data by cluster
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% cb cx cluster
% clusgsndata{1}=data.allndata(clusidx==4 | clusidx==10,:);
% clusgsndata{2}=data.allndata(clusidx==2,:);
% clusgsndata{3}=data.allndata(clusidx==6,:);
%
% clussblmean{1}=bslresp_mean(clusidx==4 | clusidx==10);
% clussblmean{2}=bslresp_mean(clusidx==2);
% clussblmean{3}=bslresp_mean(clusidx==6);
%
% clussbslresp_sd{1}=bslresp_sd(clusidx==4 | clusidx==10);
% clussbslresp_sd{2}=bslresp_sd(clusidx==2);
% clussbslresp_sd{3}=bslresp_sd(clusidx==6);
%
% clusprefdir{1}=data.allprefdir(clusidx==4 | clusidx==10,:);
% clusprefdir{2}=data.allprefdir(clusidx==2,:);
% clusprefdir{3}=data.allprefdir(clusidx==6,:);
%
% clusssds{1}=data.allssds(clusidx==4 | clusidx==10);
% clusssds{2}=data.allssds(clusidx==2);
% clusssds{3}=data.allssds(clusidx==6);

%% cDN clusters
for clusixdnum=1:length(unique(clusidx))-1
    % raster data
    clusgsndata{clusixdnum}=data.allndata(clusidx==(100+clusixdnum),:);
    % baseline
    %     clussblmean{clusixdnum}=bslresp_mean(clusidx==(100+clusixdnum));
    clussblresp{clusixdnum}=bslresps(clusidx==(100+clusixdnum),:);
    % baseline sd
    clussbslresp_sd{clusixdnum}=bslresp_sd(clusidx==(100+clusixdnum));
    % full mean
    clussfrmean{clusixdnum}=fr_mean(clusidx==(100+clusixdnum));
    % full sd
    clussfr_sd{clusixdnum}=fr_sd(clusidx==(100+clusixdnum));
    % prefered direction
    clusprefdir{clusixdnum}=data.allprefdir(clusidx==(100+clusixdnum),:);
    % ssds
    clusssds{clusixdnum}=data.allssds(clusidx==(100+clusixdnum));
    % saccade delays
    clussacRT{clusixdnum}=data.allsacdelay(clusidx==(100+clusixdnum));
    % prevalent ssds
    clusprevssd{clusixdnum}=data.allprevssd(clusidx==(100+clusixdnum));
    % ssrts
    clusmssrt{clusixdnum}=data.allmssrt_tacho(clusidx==(100+clusixdnum));
    % database info
    clusdbinfo{clusixdnum}=data.alldb(clusidx==(100+clusixdnum));
end

%% prealloc compile data
arraysz=max([sum(clusidx==101), sum(clusidx==102), sum(clusidx==103), sum(clusidx==104)]);
compgssdf=struct('clus',{'rampfallclus','sacburstclus','rampatw','fallatw'},...
    'align',struct('sac',struct('NSStrial',nan(arraysz,1301),'CStrial',nan(arraysz,1301),'NCStrial',nan(arraysz,1301),'evttimes',nan(arraysz,7)),...
    'tgt',struct('NSStrial',nan(arraysz,901),'CStrial',nan(arraysz,901),'NCStrial',nan(arraysz,901),'evttimes',nan(arraysz,7),'nnorm_nss_sdf',nan(arraysz,1301)),...
    'ssd',struct('LMCS_NSStrial',nan(arraysz,1501),'LMNCS_NSStrial',nan(arraysz,1501),'CStrial',nan(arraysz,1501),'NCStrial',nan(arraysz,1501),...
    'evttimes',nan(arraysz,7),'RPE',struct('LMCS_NSStrial',[],'LMNCS_NSStrial',[],'CStrial',[],'NCStrial',[],'PETH_REPburst_time',[],'timing',{},'amplitude',{},...
    'nxtrial_type',{},'nxtrial_pktiming',{},'nxtrial_bsline_pk',{})),...
    'corsac',struct('NSStrial',nan(arraysz,1001),'CStrial',nan(arraysz,1001),'NCStrial',nan(arraysz,1001),'evttimes',nan(arraysz,7)),...
    'rew',struct('NSStrial',nan(arraysz,1001),'CStrial',nan(arraysz,1001),'NCStrial',nan(arraysz,1001),'evttimes',nan(arraysz,7))));

trialtype={'NSStrial','CStrial','NCStrial'};
ssdtrialtype={'LMCS_NSStrial','LMNCS_NSStrial','CStrial','NCStrial'};

% data span
sac_startstop=[900 400];
tgt_startstop=[200 700];
ssd_startstop=[800 700];
corsac_startstop=[800 200];
rew_startstop=[800 200];

% colors for population plots
%     figure(1);
% close all
cc=lines(size(data.allalignmnt,1)); % one color per file
if size(cc,1)==8
    cc(8,:)=[0 0.75 0];
end

%% calculate sdf in main cluster
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clusnum=1;
for gsd=1:size(clusgsndata{clusnum},1) %compile data across all files for each condition
    %% tgt alignment ('correct_slow')
    gsdata=clusgsndata{clusnum}{gsd,2}; %this will be 3 structures containing rasters
    % and alignments tgt/tgt-CS/tgt-NCS)
    %         try
    %             gspk=clusallpk{clusnum}{gsd,2}.vis;
    %         catch nopeak
    %             gspk=0;
    %         end
    if proc_option.singlessd && iscell(clusmssrt{clusnum}{gsd,1}) &&...
            size(gsdata(1, 2).rast,1) > 5 && size(gsdata(1, 3).rast,1) > 1 &&...
            ~isempty(clusprevssd{clusnum}{gsd,1}{:})
        % Keeping NSS trials with sac latencies long enough
        % that they would have occured after a stop-signal
        if size(clussacRT{clusnum}{gsd,1}{:},2)~=size(gsdata(1, 1).rast,1)
            bugf{clusnum}(gsd,1)=clusdbinfo{clusnum}{gsd,1}.rec_id;
            %                 gsdata=[]
            %                 continue
        else
            %most prevalent SSD
            % turns out there are more empty prevssd than meet the eye:
            %                         cellfun(@(x) isempty(x{:}),clusprevssd{1, 1}(~cellfun('isempty',clusprevssd{1, 1})))
            
            [ssdbin,ssdbinval]=hist([clusssds{clusnum}{gsd,1}{1, 1};clusssds{clusnum}{gsd,1}{1, 2}]);
            ssdspread=abs(clusprevssd{clusnum}{gsd,1}{:}-max(ssdbinval(ssdbin==max(ssdbin))));
            prevssd=clusprevssd{clusnum}{gsd,1}{:}(ssdspread==min(ssdspread));
            %latency matched NSS trials
            gsdata(1, 1).rast=gsdata(1, 1).rast(clussacRT{clusnum}{gsd,1}{:}>prevssd+round(clusmssrt{clusnum}{gsd,1}{1}),:);
            
            %keep relevant SS trials (with SSD within +-3ms of prevalent SSD)
            %CS trials
            gsdata(1, 2).rast=gsdata(1, 2).rast(logical(arrayfun(@(x) sum(prevssd<=x+3 & prevssd>=x-3),...
                clusssds{clusnum}{gsd,1}{1, 1})),:);
            %NCS trials
            gsdata(1, 3).rast=gsdata(1, 3).rast(logical(arrayfun(@(x) sum(prevssd<=x+3 & prevssd>=x-3),...
                clusssds{clusnum}{gsd,1}{1, 2})),:);
        end
        %         else
        %             gsdata=[];
    end
    if ~isempty(gsdata) && clussbslresp_sd{clusnum}(gsd)~=0
        for tgtalg=1:3
            if proc_option.defaultplot && tgtalg==3
                continue
            end
            try
                rasters=gsdata(tgtalg).rast;
                if proc_option.prefdironly
                    rasters=rasters(clusprefdir{clusnum}{gsd,1}{tgtalg},:); % sorted by 3 structures containing logical indices
                end
                alignmtt=gsdata(tgtalg).alignt;
                start=alignmtt-tgt_startstop(1)-half_sixsig; stop=alignmtt+tgt_startstop(2)+half_sixsig;
                normsdf=conv_raster(rasters,conv_sigma,start,stop,clussblresp{clusnum}(gsd,:)); %normalize by baseline ,clussbslresp_sd{clusnum}(gsd)
                
                %normalize sdf by baseline activity
                %                     normsdf=(sdf-clussfrmean{clusnum}(gsd))./clussfr_sd{clusnum}(gsd);
                %                     normsdf=(sdf-clussblmean{clusnum}(gsd))./clussbslresp_sd{clusnum}(gsd);
                
                %keep ssd+ssrt time
                compgssdf(1,clusnum).align.tgt.evttimes(gsd,6:7)=round([tgt_startstop(1)+prevssd+clusmssrt{1, clusnum}{gsd, 1}{1, 1}...
                    clusmssrt{1, clusnum}{gsd, 1}{1, 3}]); %first number is ssd+ssrt, second number is tachowidth
                %keep first alignment's sdf, not normalized!
                if tgtalg==1
                    compgssdf(1,clusnum).align.tgt.nnorm_nss_sdf(gsd,:)=conv_raster(rasters,conv_sigma,start-(600-tgt_startstop(1)),stop);
                end
                %% plots (get [sdf, convrasters, convrastsem] if needed)
                %                                                  figure(1)
                %                                                  hold off
                % %                                                  patch([1:length(normsdf),fliplr(1:length(normsdf))],[normsdf-convrastsem,fliplr(normsdf+convrastsem)],'k','EdgeColor','none','FaceAlpha',0.1);
                %                                                  hold on
                %                                                  %plot sdf
                %                                                  plot(normsdf,'Color','b','LineWidth',1.8);
                %                                                  set(gca,'xtick',[1:100:1301],'xticklabel',[-200:100:700])
                %                                                  close(gcf)
                
                %% store
                compgssdf(1,clusnum).align.tgt.(trialtype{tgtalg})(gsd,:)=normsdf;
                
                %% find how early curves might diverge
                if proc_option.singlessd && tgtalg==2
                    NSStrialSDF=compgssdf(1,clusnum).align.tgt.nnorm_nss_sdf(gsd,:);
                    CStrialSDF=conv_raster(rasters,conv_sigma,start-(600-tgt_startstop(1)),stop);
                    %threshold will be 2SD above mean of differential
                    %SDF in 500ms interval before target presentation
                    BslDiff=abs(NSStrialSDF(1:600)-CStrialSDF(1:600));
                    DivThd=mean(BslDiff)+2*(std(BslDiff));
                    Diffsdfs=abs(NSStrialSDF(600-tgt_startstop(1):end)-CStrialSDF(600-tgt_startstop(1):end));
                    [DivTime,init_DivTime]=deal(tgt_startstop(1)+compgssdf(1,clusnum).align.tgt.evttimes(gsd,6));
                    if compgssdf(1,clusnum).align.tgt.evttimes(gsd,6)<tgt_startstop(2) & ... % if SSRT within range
                            logical(sum(Diffsdfs(tgt_startstop(1):end)>DivThd))    %Differential sf above threshold
                        % find max divergence
                        maxDiv=max(Diffsdfs(tgt_startstop(1):end));maxDiv_t=tgt_startstop(1)+find(Diffsdfs(tgt_startstop(1):end)==maxDiv,1);
                        [DivTime,compgssdf(1,clusnum).align.tgt.evttimes(gsd,5)]=...
                            deal(find(Diffsdfs(maxDiv_t:end)<DivThd,1,'first')-DivTime+maxDiv_t-1); %time from SSRT: end significant divergence
                        DivTime=DivTime+init_DivTime;
                        while Diffsdfs(DivTime-1)>DivThd & DivTime>1
                            DivTime=DivTime-1;
                        end
                        compgssdf(1,clusnum).align.tgt.evttimes(gsd,4)=DivTime-init_DivTime; %time from SSRT: begin significant divergence
                        % plot cluster one figures
                        if clusnum==1 & proc_option.printplots==1
                            query = ['SELECT a_file FROM sorts s INNER JOIN recordings r on s.recording_id_fk = r.recording_id WHERE sort_id = ' ...
                                num2str(clusdbinfo{clusnum}{gsd, 1}.sort_id) ];
                            recname=fetch(conn,query); recname=recname{1}(1:end-1);
                            figure('Name',[recname ' rec_id ' num2str(clusdbinfo{1}{gsd}.rec_id) ' Tgt align Single SSD'],'NumberTitle','off',...
                                'position',[2100 500 560 420]);
                            plot(NSStrialSDF(600-tgt_startstop(1):end)); hold on; plot(CStrialSDF(600-tgt_startstop(1):end));
                            axis(gca,'tight'); box off;
                            plot(DivThd*ones(1,max(get(gca,'xlim'))),'--'); plot(Diffsdfs,'-.');
                            patch([init_DivTime-2:init_DivTime+2 fliplr(init_DivTime-2:init_DivTime+2)], ...
                                reshape(repmat(get(gca,'ylim'),5,1),1,numel(get(gca,'ylim'))*5), ...
                                [0 0 0],'EdgeColor','none','FaceAlpha',0.5,'marker','d');
                            patch([tgt_startstop(1)-2:tgt_startstop(1)+2 fliplr(tgt_startstop(1)-2:tgt_startstop(1)+2)], ...
                                reshape(repmat(get(gca,'ylim'),5,1),1,numel(get(gca,'ylim'))*5), ...
                                [0 1 0],'EdgeColor','none','FaceAlpha',0.5);
                            patch([tgt_startstop(1)+prevssd-2:tgt_startstop(1)+prevssd+2 fliplr(tgt_startstop(1)+prevssd-2:tgt_startstop(1)+prevssd+2)], ...
                                reshape(repmat(get(gca,'ylim'),5,1),1,numel(get(gca,'ylim'))*5), ...
                                [1 0 0],'EdgeColor','none','FaceAlpha',0.5);
                            plot(DivTime,Diffsdfs(DivTime)+1,'*');%start
                            plot(compgssdf(1,clusnum).align.tgt.evttimes(gsd,5)+init_DivTime,Diffsdfs(compgssdf(1,clusnum).align.tgt.evttimes(gsd,5)+init_DivTime)+1,'*');%end
                            set(gca,'xtick',[1:100:1301],'xticklabel',{'-200' '-100' 'target' '100' '200' '300' '400' '500' '600' '700'});
                            title([recname ' rec_id ' num2str(clusdbinfo{1}{gsd}.rec_id) ' Tgt align Single SSD'],'Interpreter','none')
                            legend('NSS trials','CS trials','2SD > mean bsl','Differential act.')
                            % print
                            choice = questdlg('Print this one?','Print figure','Yes','No','No');
                            switch choice
                                case 'No'
                                    close(gcf);
                                case 'Yes'
                                    exportfigname=get(gcf,'Name');
                                    exportfigname=strrep(exportfigname,' ','_');
                                    savefig([exportfigname '.fig']);
                                    %print png
                                    print(gcf, '-dpng', '-noui', '-opengl','-r300', exportfigname);
                                    %print svg
                                    plot2svg([exportfigname,'.svg'],gcf, 'png');
                                    close(gcf)
                            end
                        end
                    end
                end
                
            catch norast
                compgssdf(1,clusnum).align.tgt.(trialtype{tgtalg})(gsd,:)=NaN;
            end
        end
    else
        %keep nans
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute and plot population activity (and ci) by cluster
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if proc_option.popplots
    trialtype_sdf={'NSStrial_popsdf','CStrial_popsdf','NCStrial_popsdf'};
    trialtype_ci={'NSStrial_popci','CStrial_popci','NCStrial_popci'};
    ssdtrialtype_sdf={'LMCS_NSStrial_popsdf','LMNCS_NSStrial_popsdf','CStrial_popsdf','NCStrial_popsdf'};
    ssdtrialtype_ci={'LMCS_NSStrial_popci','LMNCS_NSStrial_popci','CStrial_popci','NCStrial_popci'};
    
    CondStyles = linspecer(4);
    ClusStyles = linspecer(8);
    [sacfig,tgtfig,ssdfig,corsacfig,rewfig]=deal(nan(1,3)); %handles for figures
    
    
    %% basic plots
    if proc_option.basicplots
        %% target alignment plot
        tgtfig(clusnum)=figure('name',['Cluster' num2str(clusnum) ' target plots']);
        hold on;
        popssrt=compgssdf(1, clusnum).align.tgt.evttimes(:,6);
        %     poptacho=compgssdf(1, clusnum).align.tgt.evttimes(:,7); %if we want
        %     the curve width
        for tgtpop=1:3
            try
                %individual plots
                %                     foo=(compgssdf(1,clusnum).align.tgt.(trialtype{1}));
                %                     faa=(compgssdf(1,clusnum).align.tgt.(trialtype{2}));
                %                     for mkp=1:size(foo,1)
                %                     if ~isempty(nansum(foo(mkp,:)))
                %                       figure; plot(foo(mkp,:)); hold on;  plot(faa(mkp,:),'r'); hold off
                %                     end
                %                     end
                popsdf=nanmean(compgssdf(1,clusnum).align.tgt.(trialtype{tgtpop}));
                conv_popsdf=fullgauss_filtconv(popsdf,20,0);
                %                     popsdf(60:end-60)=conv_popsdf(60:end-60);
                [popsdf, compgssdf(1,clusnum).align.tgt.(trialtype_sdf{tgtpop})]=deal(popsdf);
                %                 [popsdf, compgssdf(1,clusnum).align.tgt.(trialtype_sdf{tgtpop})]=deal(nanmean(compgssdf(1,clusnum).align.tgt.(trialtype{tgtpop})));
                [popci, compgssdf(1,clusnum).align.tgt.(trialtype_ci{tgtpop})]=deal(nanstd(compgssdf(1,clusnum).align.tgt.(trialtype{tgtpop}))/...
                    sqrt(size(compgssdf(1,clusnum).align.tgt.(trialtype{tgtpop}),1)));
            catch
                continue
            end
            %         conv_raster(repmat(popsdf,19,1)./1000,conv_sigma);
            
            % plot SEM
            patch([1:length(popci),fliplr(1:length(popci))],...
                [popsdf-popci,fliplr(popsdf+popci)],...
                CondStyles(tgtpop,:),'EdgeColor','none','FaceAlpha',0.1);
            hold on;
            % plot SSD+SSRTs
            popssrt=popssrt(~isnan(popssrt));
            currylim=get(gca,'ylim');
            for ppssrt=1:length(popssrt)
                ssrtph=patch([popssrt(ppssrt)-1:popssrt(ppssrt)+1 fliplr(popssrt(ppssrt)-1:popssrt(ppssrt)+1)], ...
                    reshape(repmat([currylim(1) currylim(2)],3,1),1,numel(currylim)*3), ...
                    [0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.5);
            end
            lineh(tgtpop)=plot(popsdf,'LineWidth',2,'color',CondStyles(tgtpop,:));
        end
        currylim=get(gca,'ylim');
        patch([tgt_startstop(1)-2:tgt_startstop(1)+2 fliplr(tgt_startstop(1)-2:tgt_startstop(1)+2)], ...
            reshape(repmat([currylim(1) currylim(2)],5,1),1,numel(currylim)*5), ...
            [0.2 0 0.4],'EdgeColor','none','FaceAlpha',0.5);
        % pop n
        text(diff(get(gca,'xlim'))/10, currylim(2)-diff(get(gca,'ylim'))/10, ['n=' num2str(sum(clusidx==100+clusnum))])
        
        hold off;
        %% beautify plot
        set(gca,'XTick',[0:100:(tgt_startstop(2)+tgt_startstop(1))]);
        set(gca,'XTickLabel',[-tgt_startstop(1):100:tgt_startstop(2)]);
        axis(gca,'tight'); box off;
        set(gca,'Color','white','TickDir','out','FontName','Cambria','FontSize',10);
        hxlabel=xlabel(gca,'Time (ms)','FontName','Cambria','FontSize',10);
        hylabel=ylabel(gca,'Firing rate (z-score)','FontName','Cambria','FontSize',10);
        if proc_option.defaultplot
            legh=legend([lineh(1:2) ssrtph],{'No Stop Signal','Stop Signal: Cancelled','Stop Sig. Reac. Time'});
        else
            legh=legend([lineh(1:3) ssrtph],{'No Stop Signal','Stop Signal: Cancelled', 'Stop Signal: Non Cancelled','Stop Sig. Reac. Time'});
        end
        set(legh,'Interpreter','none','Location','NorthEast','Box', 'off','LineWidth',1.5,'FontName','Cambria','FontSize',9); % Interpreter prevents underscores turning character into subscript
        title(['Cluster' num2str(clusnum) ' Aligned to target'],'FontName','Cambria','FontSize',15);
    end
end

%% print plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if proc_option.printplots
    fighandles=[sacfig,tgtfig,ssdfig,corsacfig,rewfig];fighandles=fighandles(~isnan(fighandles));
    cd('E:\Data\Analysis\Countermanding\popclusters');
    
    for printfig=1:length(fighandles)
        exportfigname=get(get(get(fighandles(printfig),'CurrentAxes'),'title'),'String');
        exportfigname=strrep(exportfigname,' ','_');
        %print png
        % newpos =  get(fighandles(printfig),'Position')/60;
        % set(fighandles(printfig),'PaperUnits','inches','PaperPosition',newpos);
        print(fighandles(printfig), '-dpng', '-noui', '-opengl','-r600', exportfigname);
        % -noui stands for: suppress the printing of user interface (ui) controls.
        
        %print pdf
        %reasonably low size / good definition pdf figure (but patch transparency not supported by ghostscript to generate pdf):
        %print(fighandles(printfig), '-dpdf', '-noui', '-painters','-r600', exportfigname);
        
        %print svg
        plot2svg([exportfigname,'.svg'],fighandles(printfig), 'png'); %only vector graphic export function that preserves alpha transparency
    end
    
    %% Histograms on time of activity difference
    if proc_option.singlessd
        diffstart=compgssdf(1).align.tgt.evttimes(:,4);
        figure('Name','Timing of activity difference');
        subplot(2,1,1)
        histogram(diffstart(~isnan(diffstart)),-600:100:100); title('Start of activity difference');
        diffend=compgssdf(1).align.tgt.evttimes(:,5);
        subplot(2,1,2)
        histogram(diffend(~isnan(diffstart)),-600:100:200); title('End of activity difference');
        
        exportfigname=get(gcf,'Name');
        exportfigname=strrep(exportfigname,' ','_');
        savefig([exportfigname '.fig']);
        %print png
        print(gcf, '-dpng', '-noui', '-opengl','-r300', exportfigname);
        %print svg
        %         plot2svg([exportfigname,'.svg'],gcf, 'png');
        print(gcf, '-dsvg', '-noui', exportfigname);
        close(gcf)
    end
end
end








