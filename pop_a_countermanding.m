function pop_a_countermanding(data,recluster,conn)

%% processing options
prefdironly=0;
singlessd=1;
    if singlessd
        prefdironly=0; %otherwise we don't keep anything
    end
basicplots=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% reshape data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

allgsalignmnt=reshape({data.aligntype},size(data));
allgsprevssd=reshape({data.prevssd},size(data));allgsprevssd=allgsprevssd(:,1);
allgsssds=reshape({data.ssds},size(data));allgsssds=allgsssds(:,1);
allgssacdelay=reshape({data.sacdelay},size(data));allgssacdelay=allgssacdelay(:,1);
allgsprefdir=reshape({data.prefdiridx},size(data));
allgsndata=reshape({data.ndata},size(data));
allgsmssrt_tacho=reshape({data.allmssrt_tacho},size(data)); allgsmssrt_tacho=allgsmssrt_tacho(:,1);

% allgstasks=reshape({data.task},size(data)); allgstasks=allgstasks(:,1);
% allgspk=reshape({data.pk},size(data));
% allgs_rec_id=reshape({data.db_rec_id},size(data));
% allgsstats=reshape({data.stats},size(data));
% allgstrialidx=reshape({data.trialidx},size(data));
% allgsfname=reshape({data.fname},size(data)); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% normalize data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%convolve rasters with 200ms before saccade, 200 after saccade, 20ms kernel
%time window. Add kernel * 6 ms (see fullgauss_filtconv), e.g. 60 ms at both
% ends, which will be cut.
% allgsndata has 3 column for 3 aligntype. Each cell has 3 or 4 for different conditions

sacresps=cellfun(@(x) conv_raster(x(1,1).rast,10,x(1,1).alignt-260,x(1,1).alignt+259), allgsndata(:,1), 'UniformOutput',false); %400ms period
bslresps=cellfun(@(x) conv_raster(x(1,1).rast,10,x(1,1).alignt-720,x(1,1).alignt-1), allgsndata(:,2), 'UniformOutput',false); %600ms period

%% remove bad apples
badapl=cellfun(@(x) size(x,2)==1, sacresps);
sacresps=sacresps(~badapl,:);
sacresps=cat(1,sacresps{:});
bslresps=bslresps(~badapl,:);
bslresps=cat(1,bslresps{:});

allgsalignmnt=allgsalignmnt(~badapl,:);
allgsmssrt_tacho=allgsmssrt_tacho(~badapl,1);
%allgspk=allgspk(~badapl,:); %not needed
allgsndata=allgsndata(~badapl,:);
%allgs_rec_id=allgs_rec_id(~badapl,1); %not needed
%allgsstats=allgsstats(~badapl,1); %not needed
allgsprevssd=allgsprevssd(~badapl,:);
allgsssds=allgsssds(~badapl,:);
allgssacdelay=allgssacdelay(~badapl,:);
allgsprefdir=allgsprefdir(~badapl,:);
%allgstrialidx=allgstrialidx(~badapl,:); %not needed

%plot population
% figure;
% imagesc(sacresps);
% colormap gray
% colorbar;

%% standardize response
% z-score normalization by baseline - based on pre-target response
bslresp_mean=nanmean(bslresps');
bslresp_sd=nanstd(bslresps');
bnorm_sacresps=(sacresps-repmat(bslresp_mean',1,size(sacresps,2)))./repmat(bslresp_sd',1,size(sacresps,2));

% z-score normalization over response period (alternative method, if typically low
% baseline firing rate). Also forces clustering to operate on shapes rather than
% amplitude, by squashing response range
sacresp_mean=nanmean(sacresps');
sacresp_sd=nanstd(sacresps');
rnorm_sacresps=(sacresps-repmat(sacresp_mean',1,size(sacresps,2)))./repmat(sacresp_sd',1,size(sacresps,2));

%plot standardized population
% figure;
% imagesc(rnorm_sacresps);
% colormap gray
% colorbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Cluster population
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if recluster
    method='hclus';
    [clusidx,clustypes,clusavwf]=clus_pop(sacresps,bnorm_sacresps,rnorm_sacresps,method);
    [clusidx,sortidx]=sort(clusidx)
    
    figure
    subplot(1,2,1)
    imagesc(1:size(bnorm_sacresps,2),1:size(bnorm_sacresps,1),bnorm_sacresps)
    set(gca,'FontSize',18);
    xlabel('Time')
    ylabel('Cell Response Number')
    title('Unsorted')

    subplot(1,2,2)
    imagesc(1:size(bnorm_sacresps,2),1:size(bnorm_sacresps,1),bnorm_sacresps(sortidx,:))
    set(gca,'FontSize',18);
    xlabel('Time')
    ylabel('Cell Response Number')
    title('Sorted')
    
    % add / change unit's profile
    conn
else
    % access unit profiles
    conn
    %clusidx
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% separate data by cluster
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% cb cx cluster
% clusgsndata{1}=allgsndata(clusidx==4 | clusidx==10,:);
% clusgsndata{2}=allgsndata(clusidx==2,:);
% clusgsndata{3}=allgsndata(clusidx==6,:);
%
% clussblmean{1}=bslresp_mean(clusidx==4 | clusidx==10);
% clussblmean{2}=bslresp_mean(clusidx==2);
% clussblmean{3}=bslresp_mean(clusidx==6);
%
% clussbslresp_sd{1}=bslresp_sd(clusidx==4 | clusidx==10);
% clussbslresp_sd{2}=bslresp_sd(clusidx==2);
% clussbslresp_sd{3}=bslresp_sd(clusidx==6);
%
% clusprefdir{1}=allgsprefdir(clusidx==4 | clusidx==10,:);
% clusprefdir{2}=allgsprefdir(clusidx==2,:);
% clusprefdir{3}=allgsprefdir(clusidx==6,:);
%
% clusssds{1}=allgsssds(clusidx==4 | clusidx==10);
% clusssds{2}=allgsssds(clusidx==2);
% clusssds{3}=allgsssds(clusidx==6);

%% cDN clusters
clus1idx=clusidx==6 | clusidx==8 | clusidx==11;
clus2idx=clusidx==9;
clus3idx=clusidx==10;
% raster data
clusgsndata{1}=allgsndata(clus1idx,:);
clusgsndata{2}=allgsndata(clus2idx,:);
clusgsndata{3}=allgsndata(clus3idx,:);
% baseline
clussblmean{1}=bslresp_mean(clus1idx);
clussblmean{2}=bslresp_mean(clus2idx);
clussblmean{3}=bslresp_mean(clus3idx);
% response baseline
clussbslresp_sd{1}=bslresp_sd(clus1idx);
clussbslresp_sd{2}=bslresp_sd(clus2idx);
clussbslresp_sd{3}=bslresp_sd(clus3idx);
% prefered direction
clusprefdir{1}=allgsprefdir(clus1idx,:);
clusprefdir{2}=allgsprefdir(clus2idx,:);
clusprefdir{3}=allgsprefdir(clus3idx,:);
% ssds
clusssds{1}=allgsssds(clus1idx);
clusssds{2}=allgsssds(clus2idx);
clusssds{3}=allgsssds(clus3idx);
% saccade delays
clussacRT{1}=allgssacdelay(clus1idx);
clussacRT{2}=allgssacdelay(clus2idx);
clussacRT{3}=allgssacdelay(clus3idx);
% prevalent ssds
clusprevssd{1}=allgsprevssd(clus1idx);
clusprevssd{2}=allgsprevssd(clus2idx);
clusprevssd{3}=allgsprevssd(clus3idx);
% ssrts
clusmssrt{1}=allgsmssrt_tacho(clus1idx);
clusmssrt{2}=allgsmssrt_tacho(clus2idx);
clusmssrt{3}=allgsmssrt_tacho(clus3idx);

%% prealloc compile data
arraysz=max([sum(clus1idx), sum(clus2idx), sum(clus3idx)]);
compgssdf=struct('clus',{'rampfallclus','sacburstclus','rewrampclus'},...
    'align',struct('sac',struct('NSStrial',nan(arraysz,1301),'CStrial',nan(arraysz,1301),'NCStrial',nan(arraysz,1301),'evttimes',nan(arraysz,7)),...
                'tgt',struct('NSStrial',nan(arraysz,901),'CStrial',nan(arraysz,901),'NCStrial',nan(arraysz,901),'evttimes',nan(arraysz,7)),...
                'ssd',struct('LMCS_NSStrial',nan(arraysz,1501),'LMNCS_NSStrial',nan(arraysz,1501),'CStrial',nan(arraysz,1501),'NCStrial',nan(arraysz,1501),'evttimes',nan(arraysz,7)),...
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

%% colors for population plots
%     figure(1);
close all
cc=lines(size(allgsalignmnt,1)); % one color per file
if size(cc,1)==8
    cc(8,:)=[0 0.75 0];
end

%% calculate sdf for each outcome in each condition in each cluster (yes, that's nested loops)
% look at clusters 2,3 and 10

conv_sigma=20;
half_sixsig=conv_sigma*3; %half kernel window

for clusnum=1:3
    for gsd=1:size(clusgsndata{clusnum},1) %compile data across all files for each condition
        
        %% sac alignment ('failed_fast')
        gsdata=clusgsndata{clusnum}{gsd,1}; %this will be 3 structures containing rasters
        % and alignments sac/cancellation/ wrong sac)
        %         try
        %             gspk=clusallgspk{clusnum}{gsd,1}.sac;
        %         catch nopeak
        %             gspk=0;
        %         end
        
        % We want to compute only single SSD data -- only if there's a
        % SSRT available. Also take the opportunity to skip this alignment
        % if there aren't enough trials
        if singlessd && iscell(clusmssrt{clusnum}{gsd,1}) && size(gsdata(1, 2).rast,1) > 5 && size(gsdata(1, 3).rast,1) > 5
            % keep NC trials with  NSS trials in which a saccade would have been
            % initiated even if a stop signal had occurred, but with saccade latencies
            % greater than the stop-signal delay plus a visual-response latency.
            % We take tachomc-tachowidth/2 rather than the arbitrary 50ms
            % from Hanes et al 98 (only if tachomc>= 50)
            if size(clussacRT{clusnum}{gsd,1}{:},2)~=size(gsdata(1, 1).rast,1)
                gsdata=[]; %issue should be fixed now
            else
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
            end
        else
            gsdata=[];
        end

        if ~isempty(gsdata) && clussbslresp_sd{clusnum}(gsd)~=0
            for sacalg=1:3
                if basicplots && sacalg==2
                    continue
                end
                try
                    rasters=gsdata(sacalg).rast;
                    if prefdironly
                        rasters=rasters(clusprefdir{clusnum}{gsd,1}{sacalg},:); % sorted by 3 structures containing logical indices
                    end
                    alignmtt=gsdata(sacalg).alignt;
                    start=alignmtt-sac_startstop(1)-half_sixsig; stop=alignmtt+sac_startstop(2)+half_sixsig;
                    sdf=conv_raster(rasters,conv_sigma,start,stop);
                    
                    %normalize sdf by baseline activity
                    normsdf=(sdf-clussblmean{clusnum}(gsd))./clussbslresp_sd{clusnum}(gsd);
                    
                    %keep cue on-times
                        compgssdf(1,clusnum).align.sac.evttimes(gsd,1)=round(median(cellfun(@(x) x(1,1), gsdata(sacalg).evttime))-(alignmtt-ssd_startstop(1))); %median rew signal time
                    %keep SS time
                    if sacalg==2 || sacalg==3
                        compgssdf(1,clusnum).align.sac.evttimes(gsd,6)=compgssdf(1,clusnum).align.sac.evttimes(gsd,1)+prevssd; 
                    end
                    
                    %% plots (get [sdf, convrasters, convrastsem] if needed)
                    %                              figure
                    %                              hold off
                    %                              patch([1:length(sdf),fliplr(1:length(sdf))],[sdf-convrastsem,fliplr(sdf+convrastsem)],'k','EdgeColor','none','FaceAlpha',0.1);
                    %                              hold on
                    %                              %plot sdf
                    %                              plot(sdf,'Color','b','LineWidth',1.8);
                    %                              set(gca,'xtick',[1:100:1301],'xticklabel',[-800:100:500])
                    %                              close(gcf)
                    %
                    %% store
                    compgssdf(1,clusnum).align.sac.(trialtype{sacalg})(gsd,:)=normsdf;
                    
                catch norast
                end
            end
        else
            %keep nans
        end
        
        %% tgt alignment ('correct_slow')
        gsdata=clusgsndata{clusnum}{gsd,2}; %this will be 3 structures containing rasters
        % and alignments tgt/tgt-CS/tgt-NCS)
        %         try
        %             gspk=clusallgspk{clusnum}{gsd,2}.vis;
        %         catch nopeak
        %             gspk=0;
        %         end
        if singlessd && iscell(clusmssrt{clusnum}{gsd,1}) && size(gsdata(1, 2).rast,1) > 5 && size(gsdata(1, 3).rast,1) > 5
            % Keeping NSS trials with sac latencies long enough
            % that they would have occured after a stop-signal
                
            %most prevalent SSD
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
        else
            gsdata=[];
        end
        
        if ~isempty(gsdata) && clussbslresp_sd{clusnum}(gsd)~=0
            for tgtalg=1:3
                if basicplots && tgtalg==3
                    continue
                end
                try
                    rasters=gsdata(tgtalg).rast;
                    if prefdironly
                        rasters=rasters(clusprefdir{clusnum}{gsd,1}{tgtalg},:); % sorted by 3 structures containing logical indices
                    end
                    alignmtt=gsdata(tgtalg).alignt;
                    start=alignmtt-tgt_startstop(1)-half_sixsig; stop=alignmtt+tgt_startstop(2)+half_sixsig;
                    sdf=conv_raster(rasters,conv_sigma,start,stop);
                    
                    %normalize sdf by baseline activity
                    normsdf=(sdf-clussblmean{clusnum}(gsd))./clussbslresp_sd{clusnum}(gsd);
                    
                    %keep ssd+ssrt time
                    compgssdf(1,clusnum).align.tgt.evttimes(gsd,6:7)=round([tgt_startstop(1)+prevssd+clusmssrt{1, clusnum}{gsd, 1}{1, 1}...
                        clusmssrt{1, clusnum}{gsd, 1}{1, 3}]); %first number is ssd+ssrt, second number is tachowidth
                    
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
                catch norast
                end
            end
        else
            %keep nans
        end
        
        %% ssd alignment ('ssd')
        gsdata=clusgsndata{clusnum}{gsd,3}; %this will be 4 structures containing rasters
        % and alignment CS & NCS with corresponding
        % Lm-NSS trials (Lm-CS,Lm-NCS,CS,NCS)
        %         try
        %             gspk=clusallgspk{clusnum}{gsd,3}.vis;
        %         catch nopeak
        %             gspk=0;
        %         end
        
        if ~isempty(gsdata) && clussbslresp_sd{clusnum}(gsd)~=0
            for ssdalg=1:4
                try
                    rasters=gsdata(ssdalg).rast;
                    if prefdironly
                        rasters=rasters(clusprefdir{clusnum}{gsd,1}{ssdalg},:); % sorted by 3 structures containing logical indices
                    end
                    alignmtt=gsdata(ssdalg).alignt;
                    start=alignmtt-ssd_startstop(1)-half_sixsig; stop=alignmtt+ssd_startstop(2)+half_sixsig;
                    sdf=conv_raster(rasters,conv_sigma,start,stop);
                    
                    %normalize sdf by baseline activity
                    normsdf=(sdf-clussblmean{clusnum}(gsd))./clussbslresp_sd{clusnum}(gsd);
                    
                    %keep error time for NCS
                    if ssdalg==2 %keep reward time from latency-matched NSS trials
                        compgssdf(1,clusnum).align.ssd.evttimes(gsd,4)=round(median(cellfun(@(x) x(4,2), gsdata(ssdalg).evttime))-(alignmtt-ssd_startstop(1))); %median rew signal time
                    elseif ssdalg==4 %keep error / cue off time
                        compgssdf(1,clusnum).align.ssd.evttimes(gsd,5)=round(median(cellfun(@(x) x(5,2), gsdata(ssdalg).evttime))-(alignmtt-ssd_startstop(1))); %median error signal time
                    end
                    %% plots (get [sdf, convrasters, convrastsem] if needed)
%                     if clusnum==1 && ssdalg==4
%                              figure%(1)
%                              hold off
%                              patch([1:length(sdf),fliplr(1:length(sdf))],[sdf-convrastsem,fliplr(sdf+convrastsem)],'k','EdgeColor','none','FaceAlpha',0.1);
%                              hold on
                             %plot sdf
%                              plot(sdf,'Color','b','LineWidth',1.8);
%                              set(gca,'xtick',[1:100:1501],'xticklabel',[-800:100:700])
%                              axis(gca,'tight'); box off;
%                              close(gcf)
%                     end
                    %% store
                    compgssdf(1,clusnum).align.ssd.(ssdtrialtype{ssdalg})(gsd,:)=normsdf;
                catch norast
                end
            end
        else
            %keep nans
        end
        
        %% corrective saccade alignment ('corrsacfailed')
        gsdata=clusgsndata{clusnum}{gsd,4}; %this will be 3 structures containing rasters
        % and alignments corsac/~corsac/NCS corsac)
        %         try
        %             gspk=clusallgspk{clusnum}{gsd,4}.corsac;
        %         catch nopeak
        %             gspk=0;
        %         end
        
        if ~isempty(gsdata) && clussbslresp_sd{clusnum}(gsd)~=0
            for csacalg=1:3
                try
                    rasters=gsdata(csacalg).rast;
                    if prefdironly
                        rasters=rasters(clusprefdir{clusnum}{gsd,1}{csacalg},:); % sorted by 3 structures containing logical indices
                    end
                    alignmtt=gsdata(csacalg).alignt;
                    start=alignmtt-corsac_startstop(1)-half_sixsig; stop=alignmtt+corsac_startstop(2)+half_sixsig;
                    sdf=conv_raster(rasters,conv_sigma,start,stop);
                    
                    %normalize sdf by baseline activity
                    normsdf=(sdf-clussblmean{clusnum}(gsd))./clussbslresp_sd{clusnum}(gsd);
                    
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
                end
            end
        else
            %keep nans
        end
        
        %% reward alignment ('rewcorrect_rewslow')
        gsdata=clusgsndata{clusnum}{gsd,5}; %this will be 3 structures containing rasters
        % and alignments NSS/CS/NCS)
        %         try
        %             gspk=clusallgspk{clusnum}{gsd,5}.rew;
        %         catch nopeak
        %             gspk=0;
        %         end
        
        if ~isempty(gsdata) && clussbslresp_sd{clusnum}(gsd)~=0
            for rewalg=1:3
                try
                    rasters=gsdata(rewalg).rast;
                    if prefdironly
                        rasters=rasters(clusprefdir{clusnum}{gsd,1}{rewalg},:); % sorted by 3 structures containing logical indices
                    end
                    alignmtt=gsdata(rewalg).alignt;
                    start=alignmtt-rew_startstop(1)-half_sixsig; stop=alignmtt+rew_startstop(2)+half_sixsig;
                    sdf=conv_raster(rasters,conv_sigma,start,stop);
                    
                    %normalize sdf by baseline activity
                    normsdf=(sdf-clussblmean{clusnum}(gsd))./clussbslresp_sd{clusnum}(gsd);
                    
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

trialtype_sdf={'NSStrial_popsdf','CStrial_popsdf','NCStrial_popsdf'};
trialtype_ci={'NSStrial_popci','CStrial_popci','NCStrial_popci'};
ssdtrialtype_sdf={'LMCS_NSStrial_popsdf','LMNCS_NSStrial_popsdf','CStrial_popsdf','NCStrial_popsdf'};
ssdtrialtype_ci={'LMCS_NSStrial_popci','LMNCS_NSStrial_popci','CStrial_popci','NCStrial_popci'};

lineStyles = linspecer(4);
[sacfig,tgtfig,ssdfig,corsacfig,rewfig]=deal(nan(1,3)); %handles for figures

for clusnum=1:3
    if basicplots
    %% saccade alignment plot
    sacfig(clusnum)=figure('name',['Cluster' num2str(clusnum) ' saccade plots']);
    hold on;
    saccues=compgssdf(1, clusnum).align.sac.evttimes(:,1);
    sacsst=compgssdf(1, clusnum).align.sac.evttimes(:,6);
    for sacpop=1:3
        [popsdf, compgssdf(1,clusnum).align.sac.(trialtype_sdf{sacpop})]=deal(nanmean(compgssdf(1,clusnum).align.sac.(trialtype{sacpop})));
        [popci, compgssdf(1,clusnum).align.sac.(trialtype_ci{sacpop})]=deal(nanstd(compgssdf(1,clusnum).align.sac.(trialtype{sacpop}))/...
            sqrt(size(compgssdf(1,clusnum).align.sac.(trialtype{sacpop}),1)));
        
        patch([1:length(popci),fliplr(1:length(popci))],...
            [popsdf-popci,fliplr(popsdf+popci)],...
            lineStyles(sacpop,:),'EdgeColor','none','FaceAlpha',0.1);
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
        lineh(sacpop)=plot(popsdf,'LineWidth',2,'color',lineStyles(sacpop,:));
    end
    currylim=get(gca,'ylim');
    patch([sac_startstop(1)-2:sac_startstop(1)+2 fliplr(sac_startstop(1)-2:sac_startstop(1)+2)], ...
        reshape(repmat([currylim(1) currylim(2)],5,1),1,numel(currylim)*5), ...
        [0 0 0],'EdgeColor','none','FaceAlpha',1);
    hold off;
    %% beautify plot
    set(gca,'XTick',[0:100:(sac_startstop(2)+sac_startstop(1))]);
    set(gca,'XTickLabel',[-sac_startstop(1):100:sac_startstop(2)]);
    axis(gca,'tight'); box off;
    set(gca,'Color','white','TickDir','out','FontName','Cambria','FontSize',10);
    hxlabel=xlabel(gca,'Time (ms)','FontName','Cambria','FontSize',10);
    hylabel=ylabel(gca,'Firing rate (z-score)','FontName','Cambria','FontSize',10);
    if basicplots
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
%         popsdf=compgssdf(1,clusnum).align.tgt.(trialtype{tgtpop});
%         popsdf=conv_raster(popsdf(~isnan(nanmean(popsdf,2)),:),conv_sigma);         
%         [popsdf, compgssdf(1,clusnum).align.tgt.(trialtype_sdf{tgtpop})]=deal(popsdf);
        [popsdf, compgssdf(1,clusnum).align.tgt.(trialtype_sdf{tgtpop})]=deal(nanmean(compgssdf(1,clusnum).align.tgt.(trialtype{tgtpop})));
        [popci, compgssdf(1,clusnum).align.tgt.(trialtype_ci{tgtpop})]=deal(nanstd(compgssdf(1,clusnum).align.tgt.(trialtype{tgtpop}))/...
            sqrt(size(compgssdf(1,clusnum).align.tgt.(trialtype{tgtpop}),1)));
%         conv_raster(repmat(popsdf,19,1)./1000,conv_sigma);

       % plot SEM
        patch([1:length(popci),fliplr(1:length(popci))],...
            [popsdf-popci,fliplr(popsdf+popci)],...
            lineStyles(tgtpop,:),'EdgeColor','none','FaceAlpha',0.1);
        hold on;
        % plot SSD+SSRTs
        popssrt=popssrt(~isnan(popssrt));
        currylim=get(gca,'ylim');
        for ppssrt=1:length(popssrt)
            ssrtph=patch([popssrt(ppssrt)-1:popssrt(ppssrt)+1 fliplr(popssrt(ppssrt)-1:popssrt(ppssrt)+1)], ...
        reshape(repmat([currylim(1) currylim(2)],3,1),1,numel(currylim)*3), ...
        [0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.5);
        end
        lineh(tgtpop)=plot(popsdf,'LineWidth',2,'color',lineStyles(tgtpop,:));
    end
    currylim=get(gca,'ylim');
    patch([tgt_startstop(1)-2:tgt_startstop(1)+2 fliplr(tgt_startstop(1)-2:tgt_startstop(1)+2)], ...
        reshape(repmat([currylim(1) currylim(2)],5,1),1,numel(currylim)*5), ...
        [0.2 0 0.4],'EdgeColor','none','FaceAlpha',0.5);
    hold off;
    %% beautify plot
    set(gca,'XTick',[0:100:(tgt_startstop(2)+tgt_startstop(1))]);
    set(gca,'XTickLabel',[-tgt_startstop(1):100:tgt_startstop(2)]);
    axis(gca,'tight'); box off;
    set(gca,'Color','white','TickDir','out','FontName','Cambria','FontSize',10);
    hxlabel=xlabel(gca,'Time (ms)','FontName','Cambria','FontSize',10);
    hylabel=ylabel(gca,'Firing rate (z-score)','FontName','Cambria','FontSize',10);
    if basicplots
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
    subplot(1,2,1)
    for ssdpop=1:2:3
        % keep only files with non-null data (which had correct ssd)
        nz_idx=nansum(compgssdf(1,clusnum).align.ssd.(ssdtrialtype{ssdpop}),2)>0;
        [popsdf, compgssdf(1,clusnum).align.ssd.(ssdtrialtype_sdf{ssdpop})]=deal(nanmean(compgssdf(1,clusnum).align.ssd.(ssdtrialtype{ssdpop})(nz_idx,:)));
        [popci, compgssdf(1,clusnum).align.ssd.(ssdtrialtype_ci{ssdpop})]=deal(nanstd(compgssdf(1,clusnum).align.ssd.(ssdtrialtype{ssdpop})(nz_idx,:))/...
            sqrt(size(compgssdf(1,clusnum).align.ssd.(ssdtrialtype{ssdpop})(nz_idx,:),1)));
        
        hold on;
         if ssdpop==3 %use adequate color
            patch([1:length(popci),fliplr(1:length(popci))],...
                [popsdf-popci,fliplr(popsdf+popci)],...
                lineStyles(2,:),'EdgeColor','none','FaceAlpha',0.1);
            lineh(ssdpop)=plot(popsdf,'LineWidth',2,'color',lineStyles(2,:));
        else
            patch([1:length(popci),fliplr(1:length(popci))],...
                [popsdf-popci,fliplr(popsdf+popci)],...
                lineStyles(1,:),'EdgeColor','none','FaceAlpha',0.1);
            lineh(ssdpop)=plot(popsdf,'LineWidth',2,'color',lineStyles(1,:));
        end
    end
    currylim=get(gca,'ylim');
    patch([ssd_startstop(1)-2:ssd_startstop(1)+2 fliplr(ssd_startstop(1)-2:ssd_startstop(1)+2)], ...
        reshape(repmat([currylim(1) currylim(2)],5,1),1,numel(currylim)*5), ...
        [0.2 0 0.4],'EdgeColor','none','FaceAlpha',0.5);
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
    subplot(1,2,2)
    poperrcd=compgssdf(1, clusnum).align.ssd.evttimes(:,5);
    poprew=compgssdf(1, clusnum).align.ssd.evttimes(:,4);
    for ssdpop=2:2:4
        % keep only files with non-null data (which had correct ssd)
        nz_idx=nansum(compgssdf(1,clusnum).align.ssd.(ssdtrialtype{ssdpop}),2)>0;
        [popsdf, compgssdf(1,clusnum).align.ssd.(ssdtrialtype_sdf{ssdpop})]=deal(nanmean(compgssdf(1,clusnum).align.ssd.(ssdtrialtype{ssdpop})(nz_idx,:)));
        [popci, compgssdf(1,clusnum).align.ssd.(ssdtrialtype_ci{ssdpop})]=deal(nanstd(compgssdf(1,clusnum).align.ssd.(ssdtrialtype{ssdpop})(nz_idx,:))/...
            sqrt(size(compgssdf(1,clusnum).align.ssd.(ssdtrialtype{ssdpop})(nz_idx,:),1)));
        
        hold on;
        if ssdpop==2 %keep same color as first subplot
            patch([1:length(popci),fliplr(1:length(popci))],...
                [popsdf-popci,fliplr(popsdf+popci)],...
                lineStyles(1,:),'EdgeColor','none','FaceAlpha',0.1);
            lineh(ssdpop)=plot(popsdf,'LineWidth',2,'color',lineStyles(1,:));
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
                lineStyles(3,:),'EdgeColor','none','FaceAlpha',0.1);
            lineh(ssdpop)=plot(popsdf,'LineWidth',2,'color',lineStyles(3,:));
        end
    end
    currylim=get(gca,'ylim');
    patch([ssd_startstop(1)-2:ssd_startstop(1)+2 fliplr(ssd_startstop(1)-2:ssd_startstop(1)+2)], ...
        reshape(repmat([currylim(1) currylim(2)],5,1),1,numel(currylim)*5), ...
        [0.2 0 0.4],'EdgeColor','none','FaceAlpha',0.5);
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
    
    %% corrective saccade plots
    
    corsacfig(clusnum)=figure('name',['Cluster' num2str(clusnum) ' corrective sac plots']);
    hold on;
    
    for corsacpop=1:3
        [popsdf, compgssdf(1,clusnum).align.corsac.(trialtype_sdf{corsacpop})]=deal(nanmean(compgssdf(1,clusnum).align.corsac.(trialtype{corsacpop})));
        [popci, compgssdf(1,clusnum).align.corsac.(trialtype_ci{corsacpop})]=deal(nanstd(compgssdf(1,clusnum).align.corsac.(trialtype{corsacpop}))/...
            sqrt(size(compgssdf(1,clusnum).align.corsac.(trialtype{corsacpop}),1)));
        
        patch([1:length(popci),fliplr(1:length(popci))],...
            [popsdf-popci,fliplr(popsdf+popci)],...
            lineStyles(corsacpop,:),'EdgeColor','none','FaceAlpha',0.1);
        hold on;
        lineh(corsacpop)=plot(popsdf,'LineWidth',2,'color',lineStyles(corsacpop,:));
    end
    currylim=get(gca,'ylim');
    patch([corsac_startstop(1)-2:corsac_startstop(1)+2 fliplr(corsac_startstop(1)-2:corsac_startstop(1)+2)], ...
        reshape(repmat([currylim(1) currylim(2)],5,1),1,numel(currylim)*5), ...
        [0 0 0],'EdgeColor','none','FaceAlpha',0.5);
    hold off;
    %% beautify plot
    set(gca,'XTick',[0:100:(corsac_startstop(2)+corsac_startstop(1))]);
    set(gca,'XTickLabel',[-corsac_startstop(1):100:corsac_startstop(2)]);
    axis(gca,'tight'); box off;
    set(gca,'Color','white','TickDir','out','FontName','Cambria','FontSize',10);
    hxlabel=xlabel(gca,'Time (ms)','FontName','Cambria','FontSize',10);
    hylabel=ylabel(gca,'Firing rate (z-score)','FontName','Cambria','FontSize',10);
    legh=legend(lineh(1:3),{'No Stop Signal','Stop Signal: Cancelled', 'Stop Signal: Non Cancelled'});
    set(legh,'Interpreter','none','Location','NorthWest','Box', 'off','LineWidth',1.5,'FontName','Cambria','FontSize',9); % Interpreter prevents underscores turning character into subscript
    title(['Cluster' num2str(clusnum) ' Aligned to corrective saccade'],'FontName','Cambria','FontSize',15);
    
    
    %% reward alignment plots
    
    rewfig(clusnum)=figure('name',['Cluster' num2str(clusnum) ' reward plots']);
    hold on;
    if basicplots
        rewcdtnb=2;
    else
        rewcdtnb=3;
    end
        
    for rewpop=1:rewcdtnb
        [popsdf, compgssdf(1,clusnum).align.rew.(trialtype_sdf{rewpop})]=deal(nanmean(compgssdf(1,clusnum).align.rew.(trialtype{rewpop})));
        [popci, compgssdf(1,clusnum).align.rew.(trialtype_ci{rewpop})]=deal(nanstd(compgssdf(1,clusnum).align.rew.(trialtype{rewpop}))/...
            sqrt(size(compgssdf(1,clusnum).align.rew.(trialtype{rewpop}),1)));
        
        patch([1:length(popci),fliplr(1:length(popci))],...
            [popsdf-popci,fliplr(popsdf+popci)],...
            lineStyles(rewpop,:),'EdgeColor','none','FaceAlpha',0.1);
        hold on;
        lineh(rewpop)=plot(popsdf,'LineWidth',2,'color',lineStyles(rewpop,:));
    end
    currylim=get(gca,'ylim');
    patch([rew_startstop(1)-2:rew_startstop(1)+2 fliplr(rew_startstop(1)-2:rew_startstop(1)+2)], ...
        reshape(repmat([currylim(1) currylim(2)],5,1),1,numel(currylim)*5), ...
        [0 0 0],'EdgeColor','none','FaceAlpha',0.5);
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

%% print plots
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
    %     plot2svg([exportfigname,'.svg'],fighandles(printfig), 'png'); %only vector graphic export function that preserves alpha transparency
end


end











