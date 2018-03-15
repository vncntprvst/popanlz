function pop_a_countermanding(data,proc_option,conn)
global directory slash;
if isempty(directory)
    [directory,slash]=SetUserDir;
end

% processing options
if proc_option.singlessd
    proc_option.prefdironly=0; %otherwise we don't keep anything
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% data normalization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters
    options.sigma=10;
    options.baselineLength=500;
    options.short_wds=200;
    options.short_wde=199;
    options.long_wds=400;

% units profiles
% unitsProfile=comp_sacresp(data,options);

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

unitList=data.alldb.unit_id; %unit_ids=cellfun(@(x) x.unit_id,data.alldb);
[sorted_unit_ids,sunitid_idx]=sort(unitList);
query = ['SELECT c.profile, c.profile_type FROM clusters c WHERE cluster_id IN (' sprintf('%.0f,' ,sorted_unit_ids(1:end-1)') num2str(sorted_unit_ids(end)) ')'];
profiles = fetch(conn,query);
sunitid_revidx(sunitid_idx)=1:length(unitList);
clusidx=[profiles{sunitid_revidx,2}];
clustypes={profiles{sunitid_revidx,1}};

data.clusters=[cell2table(clustypes','VariableNames',{'Profile'})...
    array2table(clusidx','VariableNames',{'SaccadeAlignedCluster'})];

% kernel
    conv_sigma=50;
    half_sixsig=conv_sigma*3; %half kernel window
if proc_option.ssdpkalign==1
    % Instead, classify according to when peak occurs with respect to stop signal
    [clusidx,clusterIDs]=SS_reclustering(data);
else
    clusterIDs=[2 101 102 103];
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
for clusixdnum=1:length(clusterIDs)
    % raster data
    clusgsndata{clusixdnum}=data.allndata(clusidx==clusterIDs(clusixdnum),:);
    
    % baseline
    %     clussblmean{clusixdnum}=bslresp_mean(clusidx==clusterIDs(clusixdnum));
%     clussblresp{clusixdnum}=bslresps(clusidx==clusterIDs(clusixdnum),:);
    % baseline sd
%     clussbslresp_sd{clusixdnum}=bslresp_sd(clusidx==clusterIDs(clusixdnum));
    % full mean
%     clussfrmean{clusixdnum}=fr_mean(clusidx==clusterIDs(clusixdnum));
    % full sd
%     clussfr_sd{clusixdnum}=fr_sd(clusidx==clusterIDs(clusixdnum));

    clusNormFactor{clusixdnum}=data.normFactor(clusidx==clusterIDs(clusixdnum),1);

    % prefered direction
    clusprefdir{clusixdnum}=data.allprefdir(clusidx==clusterIDs(clusixdnum),:);
    % ssds
    clusssds{clusixdnum}=data.allssds(clusidx==clusterIDs(clusixdnum));
    % saccade delays
    clussacRT{clusixdnum}=data.allsacdelay(clusidx==clusterIDs(clusixdnum));
    % prevalent ssds
    clusprevssd{clusixdnum}=data.allprevssd(clusidx==clusterIDs(clusixdnum));
    % ssrts
    clusmssrt{clusixdnum}=data.allmssrt_tacho(clusidx==clusterIDs(clusixdnum));
    % database info
    clusdbinfo{clusixdnum}=data.alldb(clusidx==clusterIDs(clusixdnum));
end

%% prealloc compile data
arraysz=max([sum(clusidx==clusterIDs(1)), sum(clusidx==clusterIDs(2)), sum(clusidx==clusterIDs(3)), sum(clusidx==clusterIDs(4))]);
compgssdf=struct('clus',{'earlyFall','rampSacFall','sacBurst','allTheWay'},...%{'rampfallclus','sacburstclus','rampatw','fallatw'},...
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

%% calculate sdf for each outcome in each condition in each cluster (yes, that's nested loops)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at clusters 2,3 and 10

for clusnum=1:4
    for gsd=1:size(clusgsndata{clusnum},1) %compile data across all files for each condition
        
        %% sac alignment ('failed_fast')
        gsdata=clusgsndata{clusnum}{gsd,1}; %this will be 3 structures containing rasters
        % and alignments sac/cancellation/ wrong sac)
        %         try
        %             gspk=clusallpk{clusnum}{gsd,1}.sac;
        %         catch nopeak
        %             gspk=0;
        %         end
        
        % Here we want to compute only single SSD data -- only if there's a
        % SSRT available. Also take the opportunity to skip this alignment
        % if there aren't enough trials
        if proc_option.singlessd && iscell(clusmssrt{clusnum}{gsd,1}) && size(gsdata(1, 2).rast,1) > 1 && size(gsdata(1, 3).rast,1) > 5
            % keep NC trials with NSS trials in which a saccade would have been
            % initiated even if a stop signal had occurred, but with saccade latencies
            % greater than the stop-signal delay plus a visual-response latency.
            % We take tachomc-tachowidth/2 rather than the arbitrary 50ms
            % from Hanes et al 98 (only if tachomc>= 50)
            if size(clussacRT{clusnum}{gsd,1}{:},2)~=size(gsdata(1, 1).rast,1)
                gsdata=[]; %issue should be fixed now
            else
                try
                    %most prevalent SSD
                    [ssdbin,ssdbinval]=hist([clusssds{clusnum}{gsd,1}{1, 1};clusssds{clusnum}{gsd,1}{1, 2}]);
                    ssdspread=abs(clusprevssd{clusnum}{gsd,1}{:}-max(ssdbinval(ssdbin==max(ssdbin))));
                    prevssd=clusprevssd{clusnum}{gsd,1}{:}(ssdspread==min(ssdspread));
                    %latency matched NSS trials
                    gsdata(1, 1).rast=gsdata(1, 1).rast(clussacRT{clusnum}{gsd,1}{:}>prevssd+(max([mean(clusmssrt{clusnum}{gsd,1}{2}) 50])-clusmssrt{clusnum}{gsd,1}{3}/2) ...
                        & clussacRT{clusnum}{gsd,1}{:}<prevssd+round(clusmssrt{clusnum}{gsd,1}{1}),:);
                    
                    %keep relevant SS trials (with SSD within +-3ms of prevalent SSD)
                    %CS trials
                    gsdata(1, 2).rast=gsdata(1, 2).rast(logical(arrayfun(@(x) sum(prevssd<=x+3 & prevssd>=x-3),...
                        clusssds{clusnum}{gsd,1}{1, 1})),:);
                    %NCS trials
                    gsdata(1, 3).rast=gsdata(1, 3).rast(logical(arrayfun(@(x) sum(prevssd<=x+3 & prevssd>=x-3),...
                        clusssds{clusnum}{gsd,1}{1, 2})),:);
                catch
                    gsdata=[];
                end
            end
        else
            %still calculate most prevalent SSD
            if iscell(clusprevssd{clusnum}{gsd,1}) && ~isempty(clusprevssd{clusnum}{gsd,1}{:})
                [ssdbin,ssdbinval]=hist([clusssds{clusnum}{gsd,1}{1, 1};clusssds{clusnum}{gsd,1}{1, 2}]);
                ssdspread=abs(clusprevssd{clusnum}{gsd,1}{:}-max(ssdbinval(ssdbin==max(ssdbin))));
                prevssd=clusprevssd{clusnum}{gsd,1}{:}(ssdspread==min(ssdspread));
            else %discard?
            end
        end
        
        if ~isempty(gsdata) %&& clussbslresp_sd{clusnum}(gsd)~=0
            for sacalg=1:3
                trialtype={'NSS','CSS','NCSS'};
                if proc_option.defaultplot && sacalg==2
                    continue
                end
                try
                    rasters=gsdata(sacalg).rast;
                    if proc_option.prefdironly
                        rasters=rasters(clusprefdir{clusnum}{gsd,1}{sacalg},:); % sorted by 3 structures containing logical indices
                    end
                    alignmtt=gsdata(sacalg).alignt;
                    start=alignmtt-sac_startstop(1)-half_sixsig; stop=alignmtt+sac_startstop(2)+half_sixsig;
                    normsdf=conv_raster(rasters,conv_sigma,start,stop,clusNormFactor{clusnum}(gsd,:)); %normalize by baseline ,clussbslresp_sd{clusnum}(gsd)
                    
                    %normalize sdf by baseline activity
                    %                     normsdf=(sdf-clussfrmean{clusnum}(gsd))./clussfr_sd{clusnum}(gsd);
                    %                     normsdf=(sdf-clussblmean{clusnum}(gsd))./clussbslresp_sd{clusnum}(gsd);
                    
                    %keep cue on-times
                    compgssdf(1,clusnum).align.sac.evttimes(gsd,1)=round(median(cellfun(@(x) x(1,1), gsdata(sacalg).evttime))-(alignmtt-ssd_startstop(1))); %median rew signal time
                    %keep SS time
                    if sacalg==2 || sacalg==3
                        compgssdf(1,clusnum).align.sac.evttimes(gsd,6)=compgssdf(1,clusnum).align.sac.evttimes(gsd,1)+prevssd;
                    end
                    
                    %% plots (get [sdf, convrasters, convrastsem] if needed)
                    %                     event_times=gsdata(sacalg).evttime;
                    %                     conditions.evtsort=0; %default 0 for chronological display order
                    %                     query = ['SELECT a_file FROM sorts s INNER JOIN recordings r on s.recording_id_fk = r.recording_id WHERE sort_id = ' ...
                    %                         num2str(clusdbinfo{clusnum}{gsd, 1}.sort_id) ];
                    %                     conditions.recname=fetch(conn,query); conditions.recname=conditions.recname{1}(1:end-1);
                    %                     conditions.clus=num2str(clusnum);
                    %                     conditions.alignment='sac';
                    %                     conditions.trialtype=trialtype{sacalg};
                    %                     conditions.save=1;
                    %                     colorrasth=Raster_sdf_colorevents(rasters,event_times,start,stop,alignmtt,conv_sigma,conditions);
                    %                     drawnow;
                    %                     uiwait(colorrasth,1); %1 seconde timeout
                    %                     close(colorrasth);
                    
                    %% store
                    compgssdf(1,clusnum).align.sac.(trialtype{sacalg})(gsd,:)=normsdf;
                    
                    
                    
                    %% individual cell plots
                    if clusnum==2 && sacalg==3
                        %                                 try
                        %                                     cut_rasters=rasters(:,alignmtt-ssd_startstop(1):alignmtt+ssd_startstop(2));
                        %                                     cut_rasters(isnan(cut_rasters))=0;
                        %                                     figure;
                        %                                     subplot(2,1,1)
                        %                                     [indy, indx] = ind2sub(size(cut_rasters),find(cut_rasters));
                        %                                     plot([indx';indx'],[indy';indy'+1],'LineStyle','-'); % plot rasters
                        %                                     axis(gca,'tight'); box off;
                        %                                     subplot(2,1,2)
                        %                                     plot((normsdf))
                        %                                     ylim=get(gca,'ylim');
                        %                                     patch([ssd_startstop(1)-2:ssd_startstop(1)+2 fliplr(ssd_startstop(1)-2:ssd_startstop(1)+2)], ...
                        %                                         reshape(repmat([ylim(1) ylim(2)],5,1),1,numel(ylim)*5), ...
                        %                                         [0 0 0],'EdgeColor','none','FaceAlpha',0.5);
                        %
                        %                                     axis(gca,'tight'); box off;
                        %                                     text(diff(get(gca,'xlim'))/10, ylim(2)-diff(get(gca,'ylim'))/10, ['cell ' num2str(max(find(clusidx==(100+3),gsd)))])
                        %                                     title(['Cluster 1 - NCS - Aligned to SSD'],'FontName','Cambria','FontSize',12);
                        % %                                     exportfigname=['Clus1_NCS_SSDal_cell ' num2str(gsd)]
                        % %                                     print(gcf, '-dpng', '-noui', '-opengl','-r600', exportfigname);
                        %                                 catch
                        %                                     gsd
                        %                                 end
                    end
                    
                    
                catch norast
                    compgssdf(1,clusnum).align.sac.(trialtype{sacalg})(gsd,:)=NaN;
                end
            end
        else
            %keep nans
        end
        
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
        if ~isempty(gsdata) %&& clussbslresp_sd{clusnum}(gsd)~=0
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
                    normsdf=conv_raster(rasters,conv_sigma,start,stop,clusNormFactor{clusnum}(gsd,:)); %normalize by baseline ,clussbslresp_sd{clusnum}(gsd)
                    
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
        
        %% ssd alignment ('ssd')
        gsdata=clusgsndata{1, clusnum}{gsd, 3}; %this will be 4 structures containing rasters
        % and alignment CS & NCS with corresponding
        % Lm-NSS trials (Lm-CS,Lm-NCS,CS,NCS)
        %         try
        %             gspk=clusallpk{clusnum}{gsd,3}.vis;
        %         catch nopeak
        %             gspk=0;
        %         end
        
        if ~isempty(gsdata) %&& clussbslresp_sd{clusnum}(gsd)~=0
            for ssdalg=1:4
                try
                    rasters=gsdata(ssdalg).rast;
                    if proc_option.prefdironly
                        rasters=rasters(clusprefdir{clusnum}{gsd,1}{ssdalg},:); % sorted by 3 structures containing logical indices
                    end
                    alignmtt=gsdata(ssdalg).alignt;
                    start=alignmtt-ssd_startstop(1)-half_sixsig; stop=alignmtt+ssd_startstop(2)+half_sixsig;
                    normsdf=conv_raster(rasters,conv_sigma,start,stop,clusNormFactor{clusnum}(gsd,:)); %normalize by baseline ,clussbslresp_sd{clusnum}(gsd)
                    
                    %normalize sdf by baseline activity
                    %                     normsdf=(sdf-clussfrmean{clusnum}(gsd))./clussfr_sd{clusnum}(gsd);
                    %                     normsdf=(sdf-clussblmean{clusnum}(gsd))./clussbslresp_sd{clusnum}(gsd);
                    
                    %% plots (get [sdf, convrasters, convrastsem] if needed)
                    %                     if clusnum==1 && ssdalg==4
                    %                     event_times=gsdata(sacalg).evttime;
                    %                     conditions.evtsort=0; %default 0 for chronological display order
                    %                     query = ['SELECT a_file FROM sorts s INNER JOIN recordings r on s.recording_id_fk = r.recording_id WHERE sort_id = ' ...
                    %                         num2str(clusdbinfo{clusnum}{gsd, 1}.sort_id) ];
                    %                     conditions.recname=fetch(conn,query); conditions.recname=conditions.recname{1}(1:end-1);
                    %                     conditions.clus=num2str(clusnum);
                    %                     conditions.alignment='sac';
                    %                     conditions.save=1;
                    %                     colorrasth=Raster_sdf_colorevents(rasters,event_times,start,stop,alignmtt,conv_sigma,conditions);
                    %                     drawnow;
                    %                     uiwait(colorrasth);
                    %                     close(colorrasth);
                    %                     end
                    
                    %% store
                    compgssdf(1,clusnum).align.ssd.(ssdtrialtype{ssdalg})(gsd,:)=normsdf;
                    
                    %% individual cell plots
                    if clusnum==1 && ssdalg==4
                        %                                 try
                        %                                     cut_rasters=rasters(:,alignmtt-ssd_startstop(1):alignmtt+ssd_startstop(2));
                        %                                     cut_rasters(isnan(cut_rasters))=0;
                        %                                     figure;
                        %                                     subplot(2,1,1)
                        %                                     [indy, indx] = ind2sub(size(cut_rasters),find(cut_rasters));
                        %                                     plot([indx';indx'],[indy';indy'+1],'LineStyle','-'); % plot rasters
                        %                                     axis(gca,'tight'); box off;
                        %                                     subplot(2,1,2)
                        %                                     plot((normsdf))
                        %                                     ylim=get(gca,'ylim');
                        %                                     patch([ssd_startstop(1)-2:ssd_startstop(1)+2 fliplr(ssd_startstop(1)-2:ssd_startstop(1)+2)], ...
                        %                                         reshape(repmat([ylim(1) ylim(2)],5,1),1,numel(ylim)*5), ...
                        %                                         [0 0 0],'EdgeColor','none','FaceAlpha',0.5);
                        %
                        %                                     axis(gca,'tight'); box off;
                        %                                     text(diff(get(gca,'xlim'))/10, ylim(2)-diff(get(gca,'ylim'))/10, ['cell ' num2str(max(find(clusidx==(100+3),gsd)))])
                        %                                     title(['Cluster 1 - NCS - Aligned to SSD'],'FontName','Cambria','FontSize',12);
                        % %                                     exportfigname=['Clus1_NCS_SSDal_cell ' num2str(gsd)]
                        % %                                     print(gcf, '-dpng', '-noui', '-opengl','-r600', exportfigname);
                        %                                 catch
                        %                                     gsd
                        %                                 end
                    end
                    
                    %keep error time for NCS
                    if ssdalg==2 %keep reward time from latency-matched NSS trials
                        compgssdf(1,clusnum).align.ssd.evttimes(gsd,4)=round(median(cellfun(@(x) x(4,2), gsdata(ssdalg).evttime))-(alignmtt-ssd_startstop(1))); %median rew signal time
                    elseif ssdalg==4 %keep error / cue off time
                        compgssdf(1,clusnum).align.ssd.evttimes(gsd,5)=round(median(cellfun(@(x) x(5,2), gsdata(ssdalg).evttime))-(alignmtt-ssd_startstop(1))); %median error signal time
                    end
                    %% keep "RPE" (late burst) info
                    % number trials by trial type (some NSS trials can be in
                    % both NSS trial type)
                    if ssdalg==1 % no need to do it every time
                        compgssdf(1,clusnum).align.ssd.RPE(gsd,1).(ssdtrialtype{1})=gsdata(1, 1).trialnb;
                        compgssdf(1,clusnum).align.ssd.RPE(gsd,1).(ssdtrialtype{2})=gsdata(1, 2).trialnb;
                        compgssdf(1,clusnum).align.ssd.RPE(gsd,1).(ssdtrialtype{3})=gsdata(1, 3).trialnb;
                        compgssdf(1,clusnum).align.ssd.RPE(gsd,1).(ssdtrialtype{4})=gsdata(1, 4).trialnb;
                    end
                    
                    if ssdalg==2 | ssdalg==4
                        %                     'RPE',struct('trialtype',{},'burst',[],'timing',[],'amplitude',[],'nxtrial_type',{},'nxtrial_pktiming',[],'nxtrial_pktimediff',[],'nxtrial_bsline',[])
                        %                   % check if PETH RPE burst
                        normsdf_peaks=normsdf(ssd_startstop(1)+200:end)>std(normsdf) & [round(diff(normsdf(ssd_startstop(1)+200:end)), 2) 0]==0;
                        if sum(normsdf_peaks) %else, keep values as NaN
                            compgssdf(1,clusnum).align.ssd.RPE(gsd,1).PETH_REPburst_time(ssdalg)=...
                                find(normsdf(ssd_startstop(1)+200:end)==max(normsdf(find(normsdf_peaks)+ssd_startstop(1)+200)),1)+ssd_startstop(1)+149;
                            
                            % preallocate RPE data
                            
                            [compgssdf(1, clusnum).align.ssd.RPE(gsd, 1).timing,...
                                compgssdf(1, clusnum).align.ssd.RPE(gsd, 1).nxtrial_bsline,...
                                compgssdf(1, clusnum).align.ssd.RPE(gsd, 1).nxtrial_pktiming,...
                                compgssdf(1, clusnum).align.ssd.RPE(gsd, 1).amplitude]=deal(nan(length([gsdata.trialnb]),2));
                            
                            compgssdf(1, clusnum).align.ssd.RPE(gsd, 1).nxtrial_type=cell(length([gsdata.trialnb]),1);
                            
                            % compute trial by trial
                            for rast=1:size(rasters,1)
                                %                             close all
                                
                                convrasters=fullgauss_filtconv(rasters(rast,start:stop),conv_sigma,0).*1000;
                                %bin-sized calculation of mean and std
%                                 bsl_bins=reshape(clussblresp{clusnum}(gsd,:),conv_sigma,length(clussblresp{clusnum}(gsd,:))/conv_sigma);
%                                 meanFR=mean(nanmean(bsl_bins,2)); % should be the same as nanmean(normepochFR)
%                                 stdFR=std(nanmean(bsl_bins,2)); % better std estimate, as std(normepochFR) just overestimates std
%                                 
%                                 convrasters=(convrasters-meanFR)./stdFR;

                                  convrasters=convrasters/clusNormFactor{clusnum}(gsd);
                                    
                                % convrasters=convrasters(conv_sigma*3+1:end-3*conv_sigma);
                                
                                % figure; plot(convrasters);
                                
                                % burst or no burst
                                % peak has to be > 200ms after alignment, over PETH
                                % std, and peak (like a burst, duh)
                                convrasters_burst=convrasters(ssd_startstop(1)+150:end)>std(normsdf)+...
                                    ... %median(convrasters(ssd_startstop(1)-(gsdata(1, ssdalg).alignt-gsdata(1, ssdalg).evttime{rast}(1,1)):end)) &...
                                    min(convrasters(ssd_startstop(1)-(gsdata(1, ssdalg).alignt-gsdata(1, ssdalg).evttime{rast}(1,1)):end)) &...
                                    [round(diff(convrasters(ssd_startstop(1)+150:end)), 2) 0.1]==0;
                                %                        compgssdf(1,clusnum).align.ssd.RPE(gsd,1).
                                if logical(sum(convrasters_burst))
                                    try
                                        % burst timing
                                        % w/ respect to sac
                                        convrasters_burst_timing=find(convrasters(ssd_startstop(1)+150:end)==max(convrasters(find(convrasters_burst)+ssd_startstop(1)+150)),1)+ssd_startstop(1)+149;
                                        compgssdf(1, clusnum).align.ssd.RPE(gsd, 1).timing(gsdata(1, ssdalg).trialnb(rast),1)=...
                                            convrasters_burst_timing-(gsdata(1, ssdalg).evttime{rast}(2,1)-gsdata(1, ssdalg).alignt+ssd_startstop(1));
                                        %and w/ respect to PETH "RPE" burst
                                        
                                        compgssdf(1, clusnum).align.ssd.RPE(gsd, 1).timing(gsdata(1, ssdalg).trialnb(rast),2)=...
                                            convrasters_burst_timing-compgssdf(1,clusnum).align.ssd.RPE(gsd,1).PETH_REPburst_time(ssdalg);
                                        
                                        % burst amplitude: absolute and from median
                                        compgssdf(1, clusnum).align.ssd.RPE(gsd, 1).amplitude(gsdata(1, ssdalg).trialnb(rast),1)=...
                                            convrasters(find(convrasters(ssd_startstop(1)+150:end)==max(convrasters(find(convrasters_burst)+ssd_startstop(1)+150)),1)+ssd_startstop(1)+149);
                                        compgssdf(1, clusnum).align.ssd.RPE(gsd, 1).amplitude(gsdata(1, ssdalg).trialnb(rast),2)=...
                                            convrasters(find(convrasters(ssd_startstop(1)+150:end)==max(convrasters(find(convrasters_burst)+ssd_startstop(1)+150)),1)+ssd_startstop(1)+149)-median(convrasters);
                                        
                                        % next trial type
                                        for nxtt=1:4
                                            if ~isempty(find(compgssdf(1,clusnum).align.ssd.RPE(gsd,1).(ssdtrialtype{nxtt})==...
                                                    gsdata(1, ssdalg).trialnb(rast)+1, 1))
                                                if strcmp(ssdtrialtype{nxtt},'LMCS_NSStrial') || strcmp(ssdtrialtype{nxtt},'LMNCS_NSStrial')
                                                    compgssdf(1, clusnum).align.ssd.RPE(gsd, 1).nxtrial_type{gsdata(1, ssdalg).trialnb(rast),1}='NSStrial';
                                                else
                                                    compgssdf(1, clusnum).align.ssd.RPE(gsd, 1).nxtrial_type{gsdata(1, ssdalg).trialnb(rast),1}=ssdtrialtype{nxtt};
                                                end
                                                nxtrial_ssdalg=nxtt;
                                            end
                                        end
                                        
                                        % next trial peak timing
                                        %absolute
                                        nx_convrasters=gsdata(nxtrial_ssdalg).rast(compgssdf(1,clusnum).align.ssd.RPE(gsd,1).(ssdtrialtype{nxtrial_ssdalg})==...
                                            gsdata(1, ssdalg).trialnb(rast)+1,:);
                                        nx_convrasters=fullgauss_filtconv(nx_convrasters(start:stop),conv_sigma,0).*1000;
                                        nx_convrasters=(nx_convrasters-meanFR)./stdFR;
                                        %                                 nx_convrasters=nx_convrasters(conv_sigma*3+1:end-3*conv_sigma);
                                        
                                        %                                 figure; plot(nx_convrasters)
                                        nx_convrasters_pk=nx_convrasters(nx_convrasters(1:ssd_startstop(1)+150)>std(normsdf)+min(convrasters) &...
                                            [round(diff(nx_convrasters(1:ssd_startstop(1)+150)), 2) 0.1]==0);
                                        compgssdf(1, clusnum).align.ssd.RPE(gsd, 1).nxtrial_pktiming(gsdata(1, ssdalg).trialnb(rast),1)=...
                                            find((nx_convrasters(1:ssd_startstop(1)+150)>std(normsdf)+min(convrasters) &...
                                            [round(diff(nx_convrasters(1:ssd_startstop(1)+150)), 2) 0.1]==0),1)+find(nx_convrasters_pk==max(nx_convrasters_pk),1);
                                        
                                        % peak time diff
                                        % time diff between peak time and SS time (or presumed one for NSS)
                                        compgssdf(1, clusnum).align.ssd.RPE(gsd, 1).nxtrial_pktiming(gsdata(1, ssdalg).trialnb(rast),2)=...
                                            ssd_startstop(1)-compgssdf(1, clusnum).align.ssd.RPE(gsd, 1).nxtrial_pktiming(gsdata(1, ssdalg).trialnb(rast),1);
                                        
                                        % next trial baseline and peak amplitude
                                        cue_ontime=ssd_startstop(1)-(gsdata(nxtrial_ssdalg).alignt-...
                                            gsdata(nxtrial_ssdalg).evttime{compgssdf(1,clusnum).align.ssd.RPE(gsd,1).(ssdtrialtype{nxtrial_ssdalg})==...
                                            gsdata(1, ssdalg).trialnb(rast)+1}(1,1));
                                        compgssdf(1, clusnum).align.ssd.RPE(gsd, 1).nxtrial_bsline_pk(gsdata(1, ssdalg).trialnb(rast),1)=...
                                            nanmean(nx_convrasters(1:cue_ontime));
                                        compgssdf(1, clusnum).align.ssd.RPE(gsd, 1).nxtrial_bsline_pk(gsdata(1, ssdalg).trialnb(rast),2)=...
                                            nx_convrasters(compgssdf(1, clusnum).align.ssd.RPE(gsd, 1).nxtrial_pktiming(gsdata(1, ssdalg).trialnb(rast),1));
                                    catch
                                        %                                     rast
                                        continue
                                    end
                                end
                            end
                        end
                    end
                    
                catch norast
                    compgssdf(1,clusnum).align.ssd.(ssdtrialtype{ssdalg})(gsd,:)=NaN;
                end
            end
            %             compgssdf
        else
            %keep nans
        end
        
        %% corrective saccade alignment ('corrsacfailed')
        gsdata=clusgsndata{clusnum}{gsd,4}; %this will be 3 structures containing rasters
        % and alignments corsac/~corsac/NCS corsac)
        %         try
        %             gspk=clusallpk{clusnum}{gsd,4}.corsac;
        %         catch nopeak
        %             gspk=0;
        %         end
        
        if ~isempty(gsdata) %&& clussbslresp_sd{clusnum}(gsd)~=0
            for csacalg=1:3
                try
                    rasters=gsdata(csacalg).rast;
                    if proc_option.prefdironly
                        rasters=rasters(clusprefdir{clusnum}{gsd,1}{csacalg},:); % sorted by 3 structures containing logical indices
                    end
                    alignmtt=gsdata(csacalg).alignt;
                    start=alignmtt-corsac_startstop(1)-half_sixsig; stop=alignmtt+corsac_startstop(2)+half_sixsig;
                    normsdf=conv_raster(rasters,conv_sigma,start,stop,clusNormFactor{clusnum}(gsd,:)); %normalize by baseline ,clussbslresp_sd{clusnum}(gsd)
                    
                    %normalize sdf by baseline activity
                    %                     normsdf=(sdf-clussblmean{clusnum}(gsd))./clussbslresp_sd{clusnum}(gsd);
                    
                    %% plots (get [sdf, convrasters, convrastsem] if needed)
                    %          figure(1)
                    %          hold off
                    %          patch([1:length(sdf),fliplr(1:length(sdf))],[sdf-convrastsem,fliplr(sdf+convrastsem)],'k','EdgeColor','none','FaceAlpha',0.1);
                    %          hold on
                    %          %plot sdf
                    %          plot(sdf,'Color','b','LineWidth',1.8);
                    %          set(gca,'xtick',[1:100:1301],'xticklabel',[-800:100:500])
                    %          close(gcf)
                    
                    %% store
                    compgssdf(1,clusnum).align.corsac.(trialtype{csacalg})(gsd,:)=normsdf;
                catch norast
                    compgssdf(1,clusnum).align.corsac.(trialtype{csacalg})(gsd,:)=NaN;
                end
            end
        else
            %keep nans
        end
        
        %% reward alignment ('rewcorrect_rewslow')
        gsdata=clusgsndata{clusnum}{gsd,5}; %this will be 3 structures containing rasters
        % and alignments NSS/CS/NCS)
        %         try
        %             gspk=clusallpk{clusnum}{gsd,5}.rew;
        %         catch nopeak
        %             gspk=0;
        %         end
        
        if ~isempty(gsdata) %&& clussbslresp_sd{clusnum}(gsd)~=0
            for rewalg=1:3
                try
                    rasters=gsdata(rewalg).rast;
                    if proc_option.prefdironly
                        rasters=rasters(clusprefdir{clusnum}{gsd,1}{rewalg},:); % sorted by 3 structures containing logical indices
                    end
                    alignmtt=gsdata(rewalg).alignt;
                    start=alignmtt-rew_startstop(1)-half_sixsig; stop=alignmtt+rew_startstop(2)+half_sixsig;
                    normsdf=conv_raster(rasters,conv_sigma,start,stop,clusNormFactor{clusnum}(gsd,:)); %normalize by baseline ,clussbslresp_sd{clusnum}(gsd)
                    if max(diff(normsdf))>1
                        %                         figure; plot(normsdf);
                    end
                    %normalize sdf by baseline activity
                    %                     normsdf=(sdf-clussblmean{clusnum}(gsd))./clussbslresp_sd{clusnum}(gsd);
                    
                    %% plots (get [sdf, convrasters, convrastsem] if needed)
                    %          figure(1)
                    %          hold off
                    %          patch([1:length(sdf),fliplr(1:length(sdf))],[sdf-convrastsem,fliplr(sdf+convrastsem)],'k','EdgeColor','none','FaceAlpha',0.1);
                    %          hold on
                    %          %plot sdf
                    %          plot(sdf,'Color','b','LineWidth',1.8);
                    %          set(gca,'xtick',[1:100:1301],'xticklabel',[-800:100:500])
                    %          close(gcf)
                    
                    %% store
                    compgssdf(1,clusnum).align.rew.(trialtype{rewalg})(gsd,:)=normsdf;
                catch norast
                    compgssdf(1,clusnum).align.rew.(trialtype{rewalg})(gsd,:)=NaN;
                end
            end
        else
            %keep nans
        end
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
    
    for clusnum=1:4
        %% basic plots
        if proc_option.basicplots
            %% saccade alignment plot
            sacfig(clusnum)=figure('name',['Cluster' num2str(clusnum) ' saccade plots']);
            hold on;
            saccues=compgssdf(1, clusnum).align.sac.evttimes(:,1);
            sacsst=compgssdf(1, clusnum).align.sac.evttimes(:,6);
            for sacpop=1:3
                try
                    popsdf=nanmean(compgssdf(1,clusnum).align.sac.(trialtype{sacpop}));
                    conv_popsdf=fullgauss_filtconv(popsdf,20,0);
                    %                     popsdf(60:end-60)=conv_popsdf(60:end-60);
                    [popsdf, compgssdf(1,clusnum).align.sac.(trialtype_sdf{sacpop})]=deal(popsdf);
                    %                 [popsdf, compgssdf(1,clusnum).align.sac.(trialtype_sdf{sacpop})]=deal(nanmean(compgssdf(1,clusnum).align.sac.(trialtype{sacpop})));
                    [popci, compgssdf(1,clusnum).align.sac.(trialtype_ci{sacpop})]=deal(nanstd(compgssdf(1,clusnum).align.sac.(trialtype{sacpop}))/...
                        sqrt(size(compgssdf(1,clusnum).align.sac.(trialtype{sacpop}),1)));
                catch
                    continue
                end
                
                patch([1:length(popci),fliplr(1:length(popci))],...
                    [popsdf-popci,fliplr(popsdf+popci)],...
                    CondStyles(sacpop,:),'EdgeColor','none','FaceAlpha',0.1);
                hold on;
                % plot cues
                saccues=saccues(~isnan(saccues));
                currylim=get(gca,'ylim');
                for ppcue=1:length(saccues)
                    cueph=patch([saccues(ppcue)-1:saccues(ppcue)+1 fliplr(saccues(ppcue)-1:saccues(ppcue)+1)], ...
                        reshape(repmat([currylim(1) currylim(2)],3,1),1,numel(currylim)*3), ...
                        [0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.5);
                end
                %plot stop signals
                %         if sacpop==2 || sacpop==3
                %             sacsst=sacsst(~isnan(sacsst));
                %             currylim=get(gca,'ylim');
                %             for ppsst=1:length(sacsst)
                %                 patch([sacsst(ppsst)-1:sacsst(ppsst)+1 fliplr(sacsst(ppsst)-1:sacsst(ppsst)+1)], ...
                %             reshape(repmat([currylim(1) currylim(2)],3,1),1,numel(currylim)*3), ...
                %             lineStyles(sacpop,:),'EdgeColor','none','FaceAlpha',0.3);
                %             end
                %         end
                lineh(sacpop)=plot(popsdf,'LineWidth',2,'color',CondStyles(sacpop,:));
            end
            currylim=get(gca,'ylim');
            patch([sac_startstop(1)-2:sac_startstop(1)+2 fliplr(sac_startstop(1)-2:sac_startstop(1)+2)], ...
                reshape(repmat([currylim(1) currylim(2)],5,1),1,numel(currylim)*5), ...
                [0 0 0],'EdgeColor','none','FaceAlpha',1);
            % pop n
            text(diff(get(gca,'xlim'))/10, currylim(2)-diff(get(gca,'ylim'))/10, ['n=' num2str(sum(clusidx==100+clusnum))])
            
            hold off;
            %% beautify plot
            set(gca,'XTick',[0:100:(sac_startstop(2)+sac_startstop(1))]);
            set(gca,'XTickLabel',[-sac_startstop(1):100:sac_startstop(2)]);
            axis(gca,'tight'); box off;
            set(gca,'Color','white','TickDir','out','FontName','Cambria','FontSize',10);
            hxlabel=xlabel(gca,'Time (ms)','FontName','Cambria','FontSize',10);
            hylabel=ylabel(gca,'Firing rate (z-score)','FontName','Cambria','FontSize',10);
            
            %legend
            if proc_option.defaultplot
                legh=legend([lineh(1:2:3) cueph],{'No Stop Signal', 'Stop Signal: Non Cancelled','Target onset'});
            else
                legh=legend([lineh(1:3) cueph],{'No Stop Signal','Stop Signal: Cancelled', 'Stop Signal: Non Cancelled','Target onset'});
            end
            set(legh,'Interpreter','none','Location','SouthWest','Box','off','LineWidth',1.5,'FontName','Cambria','FontSize',9); % Interpreter prevents underscores turning character into subscript
            title(['Cluster' num2str(clusnum) ' Aligned to saccade'],'FontName','Cambria','FontSize',15);
            
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
        
        %% ssd alignment plots
        ssdfig(clusnum)=figure('name',['Cluster' num2str(clusnum) ' ssd plots']);
        % 1st plots
        subplot(2,1,1)
        for ssdpop=1:2:3
            % keep only files with non-null data (which had correct ssd)
            nz_idx=nansum(compgssdf(1,clusnum).align.ssd.(ssdtrialtype{ssdpop}),2)~=0;
            try
                popsdf=nanmean(compgssdf(1,clusnum).align.ssd.(ssdtrialtype{ssdpop})(nz_idx,:));
                conv_popsdf=fullgauss_filtconv(popsdf,20,0);
                %                 popsdf(60:end-60)=conv_popsdf(60:end-60);
                [popsdf, compgssdf(1,clusnum).align.ssd.(ssdtrialtype_sdf{ssdpop})]=deal(popsdf);
                %             [popsdf, compgssdf(1,clusnum).align.ssd.(ssdtrialtype_sdf{ssdpop})]=deal(nanmean(compgssdf(1,clusnum).align.ssd.(ssdtrialtype{ssdpop})(nz_idx,:)));
                [popci, compgssdf(1,clusnum).align.ssd.(ssdtrialtype_ci{ssdpop})]=deal(nanstd(compgssdf(1,clusnum).align.ssd.(ssdtrialtype{ssdpop})(nz_idx,:))/...
                    sqrt(size(compgssdf(1,clusnum).align.ssd.(ssdtrialtype{ssdpop})(nz_idx,:),1)));
            catch
                continue
            end
            hold on;
            if ssdpop==3 %use adequate color
                patch([1:length(popci),fliplr(1:length(popci))],...
                    [popsdf-popci,fliplr(popsdf+popci)],...
                    CondStyles(2,:),'EdgeColor','none','FaceAlpha',0.1);
                lineh(ssdpop)=plot(popsdf,'LineWidth',2,'color',CondStyles(2,:));
            else
                patch([1:length(popci),fliplr(1:length(popci))],...
                    [popsdf-popci,fliplr(popsdf+popci)],...
                    CondStyles(1,:),'EdgeColor','none','FaceAlpha',0.1);
                lineh(ssdpop)=plot(popsdf,'LineWidth',2,'color',CondStyles(1,:));
            end
        end
        currylim=get(gca,'ylim');
        patch([ssd_startstop(1)-2:ssd_startstop(1)+2 fliplr(ssd_startstop(1)-2:ssd_startstop(1)+2)], ...
            reshape(repmat([currylim(1) currylim(2)],5,1),1,numel(currylim)*5), ...
            [0.2 0 0.4],'EdgeColor','none','FaceAlpha',0.5);
        % pop n
        text(diff(get(gca,'xlim'))/10, currylim(2)-diff(get(gca,'ylim'))/10, ['n=' num2str(sum(clusidx==100+clusnum))])
        
        hold off;
        %% beautify plot
        set(gca,'XTick',[0:200:(ssd_startstop(2)+ssd_startstop(1))]);
        set(gca,'XTickLabel',[-ssd_startstop(1):200:ssd_startstop(2)]);
        axis(gca,'tight'); box off;
        set(gca,'Color','white','TickDir','out','FontName','Cambria','FontSize',10);
        hxlabel=xlabel(gca,'Time (ms)','FontName','Cambria','FontSize',10);
        hylabel=ylabel(gca,'Firing rate (z-score)','FontName','Cambria','FontSize',10);
        legh=legend(lineh(1:2:3),{'Latency Matched No Stop Signal','Stop Signal: Cancelled'});
        set(legh,'Interpreter','none','Location','NorthWest','Box', 'off','LineWidth',1.5,'FontName','Cambria','FontSize',9); % Interpreter prevents underscores turning character into subscript
        title(['Cluster' num2str(clusnum) ' NSS CS Aligned to ssd'],'FontName','Cambria','FontSize',15);
        
        %2nd plots
        subplot(2,1,2)
        poperrcd=compgssdf(1, clusnum).align.ssd.evttimes(:,5);
        poprew=compgssdf(1, clusnum).align.ssd.evttimes(:,4);
        for ssdpop=2:2:4
            % keep only files with non-null data (which had correct ssd)
            nz_idx=nansum(compgssdf(1,clusnum).align.ssd.(ssdtrialtype{ssdpop}),2)~=0;
            try
                popsdf=nanmean(compgssdf(1,clusnum).align.ssd.(ssdtrialtype{ssdpop})(nz_idx,:));
                conv_popsdf=fullgauss_filtconv(popsdf,20,0);
                %                 popsdf(60:end-60)=conv_popsdf(60:end-60);
                [popsdf, compgssdf(1,clusnum).align.ssd.(ssdtrialtype_sdf{ssdpop})]=deal(popsdf);
                %             [popsdf, compgssdf(1,clusnum).align.ssd.(ssdtrialtype_sdf{ssdpop})]=deal(nanmean(compgssdf(1,clusnum).align.ssd.(ssdtrialtype{ssdpop})(nz_idx,:)));
                [popci, compgssdf(1,clusnum).align.ssd.(ssdtrialtype_ci{ssdpop})]=deal(nanstd(compgssdf(1,clusnum).align.ssd.(ssdtrialtype{ssdpop})(nz_idx,:))/...
                    sqrt(size(compgssdf(1,clusnum).align.ssd.(ssdtrialtype{ssdpop})(nz_idx,:),1)));
            catch
                continue
            end
            hold on;
            if ssdpop==2 %keep same color as first subplot
                patch([1:length(popci),fliplr(1:length(popci))],...
                    [popsdf-popci,fliplr(popsdf+popci)],...
                    CondStyles(1,:),'EdgeColor','none','FaceAlpha',0.1);
                lineh(ssdpop)=plot(popsdf,'LineWidth',2,'color',CondStyles(1,:));
                % plot reward times
                poprew=poprew(~isnan(poprew));
                currylim=get(gca,'ylim');
                for ppssrt=1:length(poprew)
                    rewph=patch([poprew(ppssrt)-1:poprew(ppssrt)+1 fliplr(poprew(ppssrt)-1:poprew(ppssrt)+1)], ...
                        reshape(repmat([currylim(1) currylim(2)],3,1),1,numel(currylim)*3), ...
                        [0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.5);
                end
            else
                patch([1:length(popci),fliplr(1:length(popci))],...
                    [popsdf-popci,fliplr(popsdf+popci)],...
                    CondStyles(3,:),'EdgeColor','none','FaceAlpha',0.1);
                lineh(ssdpop)=plot(popsdf,'LineWidth',2,'color',CondStyles(3,:));
            end
        end
        currylim=get(gca,'ylim');
        patch([ssd_startstop(1)-2:ssd_startstop(1)+2 fliplr(ssd_startstop(1)-2:ssd_startstop(1)+2)], ...
            reshape(repmat([currylim(1) currylim(2)],5,1),1,numel(currylim)*5), ...
            [0.2 0 0.4],'EdgeColor','none','FaceAlpha',0.5);
        % pop n
        text(diff(get(gca,'xlim'))/10, currylim(2)-diff(get(gca,'ylim'))/10, ['n=' num2str(sum(clusidx==100+clusnum))])
        
        hold off;
        %% beautify plot
        set(gca,'XTick',[0:200:(ssd_startstop(2)+ssd_startstop(1))]);
        set(gca,'XTickLabel',[-ssd_startstop(1):200:ssd_startstop(2)]);
        axis(gca,'tight'); box off;
        set(gca,'Color','white','TickDir','out','FontName','Cambria','FontSize',10);
        hxlabel=xlabel(gca,'Time (ms)','FontName','Cambria','FontSize',10);
        hylabel=ylabel(gca,'Firing rate (z-score)','FontName','Cambria','FontSize',10);
        legh=legend([lineh(2:2:4) rewph],{'Lat.Matched No Stop Signal','Stop Signal: Non Cancelled','Sac trials reward'});
        set(legh,'Interpreter','none','Location','NorthWest','Box', 'off','LineWidth',1.5,'FontName','Cambria','FontSize',9); % Interpreter prevents underscores turning character into subscript
        title(['Cluster' num2str(clusnum) ' NSS NCS Aligned to ssd'],'FontName','Cambria','FontSize',15);
        
        %% control plots
        if proc_option.controlplots
            %% corrective saccade plots
            
            if clusnum==1
                corsacfig=figure('name', 'All Clusters, NCS trials, corrective sac plots'); % ['Cluster' num2str(clusnum) ' corrective sac plots']);
            else
                set(0, 'currentfigure', corsacfig);
            end
            hold on;
            
            for corsacpop=3 %1:3
                try
                    popsdf=nanmean(compgssdf(1,clusnum).align.corsac.(trialtype{corsacpop}));
                    conv_popsdf=fullgauss_filtconv(popsdf,20,0);
                    popsdf(61:end-60)=conv_popsdf;
                    [popsdf, compgssdf(1,clusnum).align.corsac.(trialtype_sdf{corsacpop})]=deal(popsdf);
                    %             [popsdf, compgssdf(1,clusnum).align.corsac.(trialtype_sdf{corsacpop})]=deal(nanmean(compgssdf(1,clusnum).align.corsac.(trialtype{corsacpop})));
                    [popci, compgssdf(1,clusnum).align.corsac.(trialtype_ci{corsacpop})]=deal(nanstd(compgssdf(1,clusnum).align.corsac.(trialtype{corsacpop}))/...
                        sqrt(size(compgssdf(1,clusnum).align.corsac.(trialtype{corsacpop}),1)));
                catch
                    continue
                end
                patch([1:length(popci),fliplr(1:length(popci))],...
                    [popsdf-popci,fliplr(popsdf+popci)],...
                    ClusStyles(clusnum+4,:),'EdgeColor','none','FaceAlpha',0.1);
                hold on;
                lineh(clusnum)=plot(popsdf,'LineWidth',2,'color',ClusStyles(clusnum+4,:));
            end
            currylim=get(gca,'ylim');
            patch([corsac_startstop(1)-2:corsac_startstop(1)+2 fliplr(corsac_startstop(1)-2:corsac_startstop(1)+2)], ...
                reshape(repmat([currylim(1) currylim(2)],5,1),1,numel(currylim)*5), ...
                [0 0 0],'EdgeColor','none','FaceAlpha',0.5);
            % pop n
            %         text(diff(get(gca,'xlim'))/10, currylim(2)-diff(get(gca,'ylim'))/10, ['n=' num2str(sum(clusidx==100+clusnum))])
            
            hold off;
            %% beautify plot
            set(gca,'XTick',[0:100:(corsac_startstop(2)+corsac_startstop(1))]);
            set(gca,'XTickLabel',[-corsac_startstop(1):100:corsac_startstop(2)]);
            axis(gca,'tight'); box off;
            set(gca,'Color','white','TickDir','out','FontName','Cambria','FontSize',10);
            hxlabel=xlabel(gca,'Time (ms)','FontName','Cambria','FontSize',10);
            hylabel=ylabel(gca,'Firing rate (z-score)','FontName','Cambria','FontSize',10);
            %     legh=legend(lineh(1:3),{'No Stop Signal','Stop Signal: Cancelled', 'Stop Signal: Non Cancelled'});
            if clusnum==4
                legh=legend(flipud(findobj(gca,'Type','line')),{'Cluster 1','Cluster 2','Cluster 3','Cluster 4'});
                set(legh,'Interpreter','none','Location','NorthWest','Box', 'off','LineWidth',1.5,'FontName','Cambria','FontSize',9); % Interpreter prevents underscores turning character into subscript
                title(['All Clusters, NCS trials, Aligned to corrective saccade'],'FontName','Cambria','FontSize',15);
                %     legend(lineh(1:3),{'No Stop Signal','Stop Signal: Cancelled', 'Stop Signal: Non Cancelled'});
                %     set(legh,'Interpreter','none','Location','NorthWest','Box', 'off','LineWidth',1.5,'FontName','Cambria','FontSize',9); % Interpreter prevents underscores turning character into subscript
                %     title(['Cluster' num2str(clusnum) ' Aligned to corrective saccade'],'FontName','Cambria','FontSize',15);
            end
            
            %% reward alignment plots
            
            rewfig(clusnum)=figure('name',['Cluster' num2str(clusnum) ' reward plots']);
            hold on;
            if proc_option.defaultplot
                rewcdtnb=2;
            else
                rewcdtnb=3;
            end
            
            for rewpop=1:rewcdtnb
                try
                    popsdf=nanmean(compgssdf(1,clusnum).align.rew.(trialtype{rewpop}));
                    conv_popsdf=fullgauss_filtconv(popsdf,20,0);
                    %                     popsdf(60:end-60)=conv_popsdf(60:end-60);
                    [popsdf, compgssdf(1,clusnum).align.rew.(trialtype_sdf{rewpop})]=deal(popsdf);
                    %             [popsdf, compgssdf(1,clusnum).align.rew.(trialtype_sdf{rewpop})]=deal(nanmean(compgssdf(1,clusnum).align.rew.(trialtype{rewpop})));
                    [popci, compgssdf(1,clusnum).align.rew.(trialtype_ci{rewpop})]=deal(nanstd(compgssdf(1,clusnum).align.rew.(trialtype{rewpop}))/...
                        sqrt(size(compgssdf(1,clusnum).align.rew.(trialtype{rewpop}),1)));
                catch
                    continue
                end
                patch([1:length(popci),fliplr(1:length(popci))],...
                    [popsdf-popci,fliplr(popsdf+popci)],...
                    CondStyles(rewpop,:),'EdgeColor','none','FaceAlpha',0.1);
                hold on;
                lineh(rewpop)=plot(popsdf,'LineWidth',2,'color',CondStyles(rewpop,:));
            end
            currylim=get(gca,'ylim');
            patch([rew_startstop(1)-2:rew_startstop(1)+2 fliplr(rew_startstop(1)-2:rew_startstop(1)+2)], ...
                reshape(repmat([currylim(1) currylim(2)],5,1),1,numel(currylim)*5), ...
                [0 0 0],'EdgeColor','none','FaceAlpha',0.5);
            % pop n
            text(diff(get(gca,'xlim'))/10, currylim(2)-diff(get(gca,'ylim'))/10, ['n=' num2str(sum(clusidx==100+clusnum))])
            
            hold off;
            %% beautify plot
            set(gca,'XTick',[0:100:(rew_startstop(2)+rew_startstop(1))]);
            set(gca,'XTickLabel',[-rew_startstop(1):100:rew_startstop(2)]);
            axis(gca,'tight'); box off;
            set(gca,'Color','white','TickDir','out','FontName','Cambria','FontSize',10);
            hxlabel=xlabel(gca,'Time (ms)','FontName','Cambria','FontSize',10);
            hylabel=ylabel(gca,'Firing rate (z-score)','FontName','Cambria','FontSize',10);
            legh=legend(lineh(1:rewcdtnb),{'No Stop Signal','Stop Signal: Cancelled','Non Cancelled - Rew Time Estimate'});
            set(legh,'Interpreter','none','Location','SouthWest','Box', 'off','LineWidth',1.5,'FontName','Cambria','FontSize',9); % Interpreter prevents underscores turning character into subscript
            title(['Cluster' num2str(clusnum) ' Aligned to reward'],'FontName','Cambria','FontSize',15);
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% trial by trial analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make two categories (use especially cluster 2): trials that see “predicted rew time” burst
% (see Schultz Reward prediction error – RPE), using a threshold estimated from PETH,
% vs no such burst (likely NSS or CS trials – maybe compared between different trial types).
% Record amplitude of said burst, and timing w/ respect to the saccade (not looking at CS trials here).
% Can adjust to delay from PETH-measured peak later in processing. See if
% a/ presence of “RPE” leads to adjustment of SS time estimate in following trial
%     (for CS trials: time diff between peak time and SS time; for NSS,
%             time diff between peak time and median SS time distribution).
%             Using median SS time estimate might blur learning-based adjustment.
%             See if one can make RT or SSRT curve over trials and calculate
%             peak time diff with that curve for a given trial.
% b/ amplitude and timing  of RPE correlate with peak timing or baseline level in following trial.
%% use ssd alignment
for clusnum=1:4
    poperrcd=compgssdf(1, clusnum).align.ssd.evttimes(:,5);
    poprew=compgssdf(1, clusnum).align.ssd.evttimes(:,4);
    
    %look at NCS trials for each unit
    
    % keep only files with non-null data (which had correct ssd)
    nz_idx=nansum(compgssdf(1,clusnum).align.ssd.(ssdtrialtype{4}),2)~=0;
    pop_ncs=compgssdf(1,clusnum).align.ssd.(ssdtrialtype{4})(nz_idx,:);
    try
        popsdf=nanmean(compgssdf(1,clusnum).align.ssd.(ssdtrialtype{4})(nz_idx,:));
        conv_popsdf=fullgauss_filtconv(popsdf,20,0);
        %         popsdf(60:end-60)=conv_popsdf(60:end-60);
    catch
        continue
    end
    %find trials with RPE rebound above threshold
    
end
%% beautify plot
%     set(gca,'XTick',[0:200:(ssd_startstop(2)+ssd_startstop(1))]);
%     set(gca,'XTickLabel',[-ssd_startstop(1):200:ssd_startstop(2)]);
%     axis(gca,'tight'); box off;
%     set(gca,'Color','white','TickDir','out','FontName','Cambria','FontSize',10);
%     hxlabel=xlabel(gca,'Time (ms)','FontName','Cambria','FontSize',10);
%     hylabel=ylabel(gca,'Firing rate (z-score)','FontName','Cambria','FontSize',10);
%     legh=legend([lineh(2:2:4) rewph],{'Lat.Matched No Stop Signal','Stop Signal: Non Cancelled','Sac trials reward'});
%     set(legh,'Interpreter','none','Location','NorthWest','Box', 'off','LineWidth',1.5,'FontName','Cambria','FontSize',9); % Interpreter prevents underscores turning character into subscript
%     title(['Cluster' num2str(clusnum) ' NSS NCS Aligned to ssd'],'FontName','Cambria','FontSize',15);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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








