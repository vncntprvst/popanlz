function pop_a_countermanding(data,recluster,conn)

%% processing options
prefdironly=0;
singlessd=1;
if singlessd
    prefdironly=0; %otherwise we don't keep anything
end
popplots=1;
printplots=0;
defaultplot=1;
basicplots=1;
controlplots=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% normalize data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1/ convolve rasters with 200ms before saccade, 200 after saccade, 20ms kernel
%time window. Add kernel * 6 ms (see fullgauss_filtconv), e.g. 60 ms at both
% ends, which will be cut.
% data.allgsndata has 3 column for 3 aligntype. Each cell has 3 or 4 for different conditions
sigma=10;
sacresps=cellfun(@(x) conv_raster(x(1,1).rast,sigma,x(1,1).alignt-(200+sigma*3),x(1,1).alignt+(199+sigma*3)), data.allgsndata(:,1), 'UniformOutput',false); %400ms period
bslresps=cellfun(@(x) conv_raster(x(1,1).rast,sigma,x(1,1).alignt-(500+sigma*3),x(1,1).alignt+(sigma*3-1)), data.allgsndata(:,2), 'UniformOutput',false); %500ms period
fullresps=cellfun(@(x) conv_raster(x(1,1).rast,sigma,1,size(x(1,1).rast,2)), data.allgsndata(:,2), 'UniformOutput',false); %full response

%% remove bad apples
badapl=cellfun(@(x) size(x,2)==1, sacresps);
sacresps=sacresps(~badapl,:);
sacresps=cat(1,sacresps{:});
bslresps=bslresps(~badapl,:);
bslresps=cat(1,bslresps{:});
fullresps=fullresps(~badapl,:);
% fullresps=cat(1,fullresps{:});

data.allgsalignmnt=data.allgsalignmnt(~badapl,:);
data.allgsmssrt_tacho=data.allgsmssrt_tacho(~badapl,1);
%allgspk=allgspk(~badapl,:); %not needed
data.allgsndata=data.allgsndata(~badapl,:);
%allgs_rec_id=allgs_rec_id(~badapl,1); %not needed
%allgsstats=allgsstats(~badapl,1); %not needed
data.allgsprevssd=data.allgsprevssd(~badapl,:);
data.allgsssds=data.allgsssds(~badapl,:);
data.allgssacdelay=data.allgssacdelay(~badapl,:);
data.allgsprefdir=data.allgsprefdir(~badapl,:);
%allgstrialidx=allgstrialidx(~badapl,:); %not needed
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
%     align=data.allgsndata{top_drop(topfig), 3}(4).alignt;
%     rasters=((data.allgsndata{top_drop(topfig), 3}(4).rast(:,align-800:align+800)));
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cluster population
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
unit_ids=cellfun(@(x) x.unit_id,data.alldb);
if recluster
    method='hclus';
    [clusidx,clustypes,clusavwf]=clus_pop(sacresps,bnorm_sacresps,rnorm_sacresps,method);
    [~,sortidx]=sort(clusidx);
    clusid=unique(clusidx);
    
    %population raster
    poprasthm=figure('name','population raster')
    subplot(1,20,1:9)
    imagesc(1:size(bnorm_sacresps,2),1:size(bnorm_sacresps,1),bnorm_sacresps)
    set(gca,'FontSize',18);
    xlabel('Time')
    ylabel('Neuron #')
    title('Unsorted')
    subplot(1,20,12:19)
    imagesc(1:size(bnorm_sacresps,2),1:size(bnorm_sacresps,1),bnorm_sacresps(sortidx,:))
    set(gca,'FontSize',18);
    xlabel('Time')
    ylabel('Neuron #')
    title('Sorted by cluster')
    rangesph=subplot(1,20,20)
    clusrange=[zeros(sum(clusidx==clusid(1)),5);ones(sum(clusidx==clusid(2)),5);...
        zeros(sum(clusidx==clusid(3)),5);ones(sum(clusidx==clusid(4)),5);zeros(sum(clusidx==clusid(5)),5)];
    imagesc(clusrange)
    set(gca,'XTick', [],'XTickLabel',[],'YTick', [],'YTickLabel',[]);
    % cd('E:\BoxSync\Box Sync\Home Folder vp35\Sync\SommerLab\projects\countermanding\popclusters')
    % exportfigname='population raster';
    % %     print(gcf, '-dpng', '-noui', '-opengl','-r600', exportfigname);
    % plot2svg([exportfigname,'.svg'],gcf, 'png');
    %
    % clusters mean response
    figure('name','clusters mean response')
    cmrtitles={'Unsorted','ramp & fall','burst','ramp all the way','ramp down'};
    for mclussp=1:length(unique(clusidx))
        subplot(5,1,mclussp)
        plot(clusavwf(mclussp,:));
        xlabel('Time')
        ylabel('Norm. Firing rate')
        title(cmrtitles{mclussp})
        set(gca,'xtick',1:100:clusavwf(mclussp,:),'xticklabel',1:100:clusavwf(mclussp,:),'TickDir','out');
        set(gca,'Color','white','FontSize',18,'FontName','calibri');
        axis(gca,'tight'); box off;
    end
    
    % add / change unit's profile
    success = addProfile(clustypes, clusidx, unit_ids, conn )
    
else
    % access unit profiles
    [sorted_unit_ids,sunitid_idx]=sort(unit_ids);
    query = ['SELECT c.profile, c.profile_type FROM clusters c WHERE cluster_id IN (' sprintf('%.0f,' ,sorted_unit_ids(1:end-1)') num2str(sorted_unit_ids(end)) ')'];
    profiles = fetch(conn,query);
    sunitid_revidx(sunitid_idx)=1:length(unit_ids);
    clusidx=[profiles{sunitid_revidx,2}];
    clustypes={profiles{sunitid_revidx,1}};
end


%% Instead, classify according to when peak occurs with respect to stop signal
% data span
ssd_startstop=[800 700];

% kernel
conv_sigma=50;
half_sixsig=conv_sigma*3; %half kernel window

for gsd=1:size(data.allgsndata,1)
    if clusidx(gsd)~=-1
        try
        rasters=data.allgsndata{gsd,3}(3).rast;
        alignmtt=data.allgsndata{gsd,3}(3).alignt;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Separate data by cluster
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% cb cx cluster
% clusgsndata{1}=data.allgsndata(clusidx==4 | clusidx==10,:);
% clusgsndata{2}=data.allgsndata(clusidx==2,:);
% clusgsndata{3}=data.allgsndata(clusidx==6,:);
%
% clussblmean{1}=bslresp_mean(clusidx==4 | clusidx==10);
% clussblmean{2}=bslresp_mean(clusidx==2);
% clussblmean{3}=bslresp_mean(clusidx==6);
%
% clussbslresp_sd{1}=bslresp_sd(clusidx==4 | clusidx==10);
% clussbslresp_sd{2}=bslresp_sd(clusidx==2);
% clussbslresp_sd{3}=bslresp_sd(clusidx==6);
%
% clusprefdir{1}=data.allgsprefdir(clusidx==4 | clusidx==10,:);
% clusprefdir{2}=data.allgsprefdir(clusidx==2,:);
% clusprefdir{3}=data.allgsprefdir(clusidx==6,:);
%
% clusssds{1}=data.allgsssds(clusidx==4 | clusidx==10);
% clusssds{2}=data.allgsssds(clusidx==2);
% clusssds{3}=data.allgsssds(clusidx==6);

%% cDN clusters
for clusixdnum=1:length(unique(clusidx))-1
    % raster data
    clusgsndata{clusixdnum}=data.allgsndata(clusidx==(100+clusixdnum),:);
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
    clusprefdir{clusixdnum}=data.allgsprefdir(clusidx==(100+clusixdnum),:);
    % ssds
    clusssds{clusixdnum}=data.allgsssds(clusidx==(100+clusixdnum));
    % saccade delays
    clussacRT{clusixdnum}=data.allgssacdelay(clusidx==(100+clusixdnum));
    % prevalent ssds
    clusprevssd{clusixdnum}=data.allgsprevssd(clusidx==(100+clusixdnum));
    % ssrts
    clusmssrt{clusixdnum}=data.allgsmssrt_tacho(clusidx==(100+clusixdnum));
end

%% prealloc compile data
arraysz=max([sum(clusidx==101), sum(clusidx==102), sum(clusidx==103), sum(clusidx==104)]);
compgssdf=struct('clus',{'rampfallclus','sacburstclus','rampatw','fallatw'},...
    'align',struct('sac',struct('NSStrial',nan(arraysz,1301),'CStrial',nan(arraysz,1301),'NCStrial',nan(arraysz,1301),'evttimes',nan(arraysz,7)),...
    'tgt',struct('NSStrial',nan(arraysz,901),'CStrial',nan(arraysz,901),'NCStrial',nan(arraysz,901),'evttimes',nan(arraysz,7)),...
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

%% colors for population plots
%     figure(1);
% close all
cc=lines(size(data.allgsalignmnt,1)); % one color per file
if size(cc,1)==8
    cc(8,:)=[0 0.75 0];
end

%% calculate sdf for each outcome in each condition in each cluster (yes, that's nested loops)
% look at clusters 2,3 and 10



for clusnum=1:4
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
            gsdata=[];
        end
        
        if ~isempty(gsdata) && clussbslresp_sd{clusnum}(gsd)~=0
            for sacalg=1:3
                if defaultplot && sacalg==2
                    continue
                end
                try
                    rasters=gsdata(sacalg).rast;
                    if prefdironly
                        rasters=rasters(clusprefdir{clusnum}{gsd,1}{sacalg},:); % sorted by 3 structures containing logical indices
                    end
                    alignmtt=gsdata(sacalg).alignt;
                    start=alignmtt-sac_startstop(1)-half_sixsig; stop=alignmtt+sac_startstop(2)+half_sixsig;
                    normsdf=conv_raster(rasters,conv_sigma,start,stop,clussblresp{clusnum}(gsd,:)); %normalize by baseline ,clussbslresp_sd{clusnum}(gsd)
                    
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
                if defaultplot && tgtalg==3
                    continue
                end
                try
                    rasters=gsdata(tgtalg).rast;
                    if prefdironly
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
                    normsdf=conv_raster(rasters,conv_sigma,start,stop,clussblresp{clusnum}(gsd,:)); %normalize by baseline ,clussbslresp_sd{clusnum}(gsd)
                    
                    %normalize sdf by baseline activity
                    %                     normsdf=(sdf-clussfrmean{clusnum}(gsd))./clussfr_sd{clusnum}(gsd);
                    %                     normsdf=(sdf-clussblmean{clusnum}(gsd))./clussbslresp_sd{clusnum}(gsd);
                    
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
                    
                            %% individual cell plots
                            if clusnum==3 && ssdalg==4
                                try
                                    cut_rasters=rasters(:,alignmtt-ssd_startstop(1):alignmtt+ssd_startstop(2));
                                    figure;
                                    subplot(2,1,1)
                                    [indy, indx] = ind2sub(size(cut_rasters),find(cut_rasters));
                                    plot([indx';indx'],[indy';indy'+1],'LineStyle','-'); % plot rasters
                                    axis(gca,'tight'); box off;
                                    subplot(2,1,2)
                                    plot((normsdf))
                                    ylim=get(gca,'ylim');
                                    patch([ssd_startstop(1)-2:ssd_startstop(1)+2 fliplr(ssd_startstop(1)-2:ssd_startstop(1)+2)], ...
                                        reshape(repmat([ylim(1) ylim(2)],5,1),1,numel(ylim)*5), ...
                                        [0 0 0],'EdgeColor','none','FaceAlpha',0.5);
                                    
                                    axis(gca,'tight'); box off;
                                    text(diff(get(gca,'xlim'))/10, ylim(2)-diff(get(gca,'ylim'))/10, ['cell ' num2str(gsd)])
                                    title(['Cluster 1 - NCS - Aligned to SSD'],'FontName','Cambria','FontSize',12);
                                    exportfigname=['Clus1_NCS_SSDal_cell ' num2str(gsd)]
                                    print(gcf, '-dpng', '-noui', '-opengl','-r600', exportfigname);
                                catch
                                    gsd
                                end
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
                            bsl_bins=reshape(clussblresp{clusnum}(gsd,:),conv_sigma,length(clussblresp{clusnum}(gsd,:))/conv_sigma);
                            meanFR=mean(nanmean(bsl_bins,2)); % should be the same as nanmean(normepochFR)
                            stdFR=std(nanmean(bsl_bins,2)); % better std estimate, as std(normepochFR) just overestimates std

                            convrasters=(convrasters-meanFR)./stdFR;
%                             convrasters=convrasters(conv_sigma*3+1:end-3*conv_sigma);
                            
%                             figure; plot(convrasters);
                            
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
                    normsdf=conv_raster(rasters,conv_sigma,start,stop,clussblresp{clusnum}(gsd,:)); %normalize by baseline ,clussbslresp_sd{clusnum}(gsd)
                    
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
                    normsdf=conv_raster(rasters,conv_sigma,start,stop,clussblresp{clusnum}(gsd,:)); %normalize by baseline ,clussbslresp_sd{clusnum}(gsd)
                    
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
if popplots
    trialtype_sdf={'NSStrial_popsdf','CStrial_popsdf','NCStrial_popsdf'};
    trialtype_ci={'NSStrial_popci','CStrial_popci','NCStrial_popci'};
    ssdtrialtype_sdf={'LMCS_NSStrial_popsdf','LMNCS_NSStrial_popsdf','CStrial_popsdf','NCStrial_popsdf'};
    ssdtrialtype_ci={'LMCS_NSStrial_popci','LMNCS_NSStrial_popci','CStrial_popci','NCStrial_popci'};
    
    lineStyles = linspecer(4);
    [sacfig,tgtfig,ssdfig,corsacfig,rewfig]=deal(nan(1,3)); %handles for figures
    
    for clusnum=1:4
        if basicplots
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
            if defaultplot
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
            if defaultplot
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
        subplot(1,2,2)
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

        if controlplots
            %% corrective saccade plots
            
%             if clusnum==1
%                 corsacfig=figure('name', 'All Clusters, NCS trials, corrective sac plots'); % ['Cluster' num2str(clusnum) ' corrective sac plots']);
%             else
%                 set(0, 'currentfigure', corsacfig);
%             end
%             hold on;
%             
%             for corsacpop=3 %1:3
%                 try
%                     popsdf=nanmean(compgssdf(1,clusnum).align.corsac.(trialtype{corsacpop}));
%                     conv_popsdf=fullgauss_filtconv(popsdf,20,0);
%                     popsdf(60:end-60)=conv_popsdf(60:end-60);
%                     [popsdf, compgssdf(1,clusnum).align.corsac.(trialtype_sdf{corsacpop})]=deal(popsdf);
%                     %             [popsdf, compgssdf(1,clusnum).align.corsac.(trialtype_sdf{corsacpop})]=deal(nanmean(compgssdf(1,clusnum).align.corsac.(trialtype{corsacpop})));
%                     [popci, compgssdf(1,clusnum).align.corsac.(trialtype_ci{corsacpop})]=deal(nanstd(compgssdf(1,clusnum).align.corsac.(trialtype{corsacpop}))/...
%                         sqrt(size(compgssdf(1,clusnum).align.corsac.(trialtype{corsacpop}),1)));
%                 catch
%                     continue
%                 end
%                 patch([1:length(popci),fliplr(1:length(popci))],...
%                     [popsdf-popci,fliplr(popsdf+popci)],...
%                     lineStyles(clusnum+4,:),'EdgeColor','none','FaceAlpha',0.1);
%                 hold on;
%                 lineh(clusnum)=plot(popsdf,'LineWidth',2,'color',lineStyles(clusnum+4,:));
%             end
%             currylim=get(gca,'ylim');
%             patch([corsac_startstop(1)-2:corsac_startstop(1)+2 fliplr(corsac_startstop(1)-2:corsac_startstop(1)+2)], ...
%                 reshape(repmat([currylim(1) currylim(2)],5,1),1,numel(currylim)*5), ...
%                 [0 0 0],'EdgeColor','none','FaceAlpha',0.5);
%             % pop n
%             %         text(diff(get(gca,'xlim'))/10, currylim(2)-diff(get(gca,'ylim'))/10, ['n=' num2str(sum(clusidx==100+clusnum))])
%             
%             hold off;
%             %% beautify plot
%             set(gca,'XTick',[0:100:(corsac_startstop(2)+corsac_startstop(1))]);
%             set(gca,'XTickLabel',[-corsac_startstop(1):100:corsac_startstop(2)]);
%             axis(gca,'tight'); box off;
%             set(gca,'Color','white','TickDir','out','FontName','Cambria','FontSize',10);
%             hxlabel=xlabel(gca,'Time (ms)','FontName','Cambria','FontSize',10);
%             hylabel=ylabel(gca,'Firing rate (z-score)','FontName','Cambria','FontSize',10);
%             %     legh=legend(lineh(1:3),{'No Stop Signal','Stop Signal: Cancelled', 'Stop Signal: Non Cancelled'});
%             if clusnum==4
%                 legh=legend(flipud(findobj(gca,'Type','line')),{'Cluster 1','Cluster 2','Cluster 3','Cluster 4'});
%                 set(legh,'Interpreter','none','Location','NorthWest','Box', 'off','LineWidth',1.5,'FontName','Cambria','FontSize',9); % Interpreter prevents underscores turning character into subscript
%                 title(['All Clusters, NCS trials, Aligned to corrective saccade'],'FontName','Cambria','FontSize',15);
%                 %     legend(lineh(1:3),{'No Stop Signal','Stop Signal: Cancelled', 'Stop Signal: Non Cancelled'});
%                 %     set(legh,'Interpreter','none','Location','NorthWest','Box', 'off','LineWidth',1.5,'FontName','Cambria','FontSize',9); % Interpreter prevents underscores turning character into subscript
%                 %     title(['Cluster' num2str(clusnum) ' Aligned to corrective saccade'],'FontName','Cambria','FontSize',15);
%             end
            
            %% reward alignment plots
            
            rewfig(clusnum)=figure('name',['Cluster' num2str(clusnum) ' reward plots']);
            hold on;
            if defaultplot
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
                    lineStyles(rewpop,:),'EdgeColor','none','FaceAlpha',0.1);
                hold on;
                lineh(rewpop)=plot(popsdf,'LineWidth',2,'color',lineStyles(rewpop,:));
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
if printplots
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
end
end











