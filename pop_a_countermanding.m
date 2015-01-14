function pop_a_countermanding(allgsalignmnt,allgsmssrt,allgspk,allgsndata,...
    allgs_rec_id,allgsstats,allgsprevssd,allgsssds,allgssacdelay,allgsprefdir)

%% processing options
prefdironly=1;
singlessd=1;

%% cluster analysis of population
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
allgsmssrt=allgsmssrt(~badapl,1);
allgspk=allgspk(~badapl,:);
allgsndata=allgsndata(~badapl,:);
allgs_rec_id=allgs_rec_id(~badapl,1);
allgsstats=allgsstats(~badapl,1);
allgsprevssd=allgsprevssd(~badapl,:);
allgsssds=allgsssds(~badapl,:);
allgssacdelay=allgssacdelay(~badapl,:);
allgsprefdir=allgsprefdir(~badapl,:);

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

% %% k-means clustering
% % define seeds
% seeds=cellfun(@(x) mean(x(1,100:200))-mean(x(1,200:300)), mat2cell(bnorm_sacresps,ones(size(bnorm_sacresps,1),1)));
% wavedropseed=seeds==max(seeds);
% waveburstseed=seeds==min(seeds);
% waveflatseed=seeds==min(abs(seeds));
% seeds=[bnorm_sacresps(wavedropseed,:);...
%     bnorm_sacresps(waveburstseed,:);...
%     bnorm_sacresps(waveflatseed,:)];
%
% %% plot seeds
% figure
% plot(seeds(1,:))
% hold on
% plot(seeds(2,:),'r')
% plot(seeds(3,:),'g')
%
% %% calculate k-means
% % at random
% % [IDX,C,sumd,D]=kmeans(sacresps,3,'dist','city','display','iter');
% % seeded
% [kMidx,kMeansClus,sumd,D]=kmeans(bnorm_sacresps,3,'dist','city','start',seeds,'display','iter');
%
% %% plot means
% figure
% plot(kMeansClus(1,:),'b')
% hold on
% plot(kMeansClus(2,:),'r')
% plot(kMeansClus(3,:),'g')
%
% %% k-means scatterplot
% figure
% scatter3(D(kMidx==1,1),D(kMidx==1,2),D(kMidx==1,3),40,'b.')
% hold on
% scatter3(D(kMidx==2,1),D(kMidx==2,2),D(kMidx==2,3),40,'r.')
% scatter3(D(kMidx==3,1),D(kMidx==3,2),D(kMidx==3,3),40,'g.')
%
% % Fit Gaussian Mixture Model using the k-means centers as the initial conditions
% % We only have mean initial conditions from the k-means algorithm, so we
% % can specify some arbitrary initial variance and mixture weights.
% % gmInitialVariance = 0.1;
% % initialSigma = cat(3,gmInitialVariance,gmInitialVariance,gmInitialVariance);
% % % Initial weights are set at 50%
% % initialWeights = [0.5 0.5 0.5];
% % % Initial condition structure for the gmdistribution.fit function
% % S.mu = kMeansClus;
% % S.Sigma = initialSigma;
% % S.PComponents = initialWeights;
%
% %% plot sdfs from each k-means cluster
% clus1=find(kMidx==1); clus2=find(kMidx==2); clus3=find(kMidx==3);
% figure('name','cluster1')
% subplotdim=[ceil(length(clus1)/2)-(2*floor(length(clus1)/10)),2+floor(length(clus1)/10)];
% for sacplot=1:length(clus1)
%     subplot(subplotdim(1),subplotdim(2),sacplot)
%     plot(bnorm_sacresps(clus1(sacplot),:));
%     ylim=get(gca,'ylim');
%     set(gca,'ylim',[0 ylim(2)]);
%     text(20,10,['sacplot ' num2str(clus1(sacplot))]);
% end
% figure('name','cluster2')
% subplotdim=[ceil(length(clus2)/2)-(2*floor(length(clus2)/10)),2+floor(length(clus2)/10)];
% for sacplot=1:length(clus2)
%     subplot(subplotdim(1),subplotdim(2),sacplot)
%     plot(bnorm_sacresps(clus2(sacplot),:));
%     ylim=get(gca,'ylim');
%     set(gca,'ylim',[0 ylim(2)]);
%     text(20,10,['sacplot ' num2str(clus2(sacplot))]);
% end
% figure('name','cluster3')
% subplotdim=[ceil(length(clus3)/2)-(2*floor(length(clus3)/10)),2+floor(length(clus3)/10)];
% for sacplot=1:length(clus3)
%     subplot(subplotdim(1),subplotdim(2),sacplot)
%     plot(bnorm_sacresps(clus3(sacplot),:));
%     ylim=get(gca,'ylim');
%     set(gca,'ylim',[0 ylim(2)]);
%     text(20,10,['sacplot ' num2str(clus3(sacplot))]);
% end
%
% %% PCA clustering
% [coeffs,PrComps,latent] = pca(bnorm_sacresps);
% % D=coeffs(:,1:8)';
% % figure
% % plot(coeffs(:,1),'.','MarkerSize',.5)
% % hold on
% % plot(coeffs(:,2),'.','MarkerSize',.5)
% % plot(coeffs(:,3),'.','MarkerSize',.5)
% % plot(coeffs(:,4),'.','MarkerSize',.5)
%
% %% scatterplot
% % in 2D
% figure
% scatter(PrComps(:,1), PrComps(:,2), 'k.');
% xlabel('PC 1'); ylabel('PC 2')
% % hold on
% % scatter(PrComps([6,11,21,24,29,32],1),PrComps([6,11,21,24,29,32],2),80,'r','filled')
% % scatter(PrComps([7,14,17,25,27],1),PrComps([7,14,17,25,27],2),80,'g','filled')
%
% %% isolate clusters
% FirstPrComps=[PrComps(:,1),PrComps(:,2)];
% try
%     gmm_fit_clusters = gmdistribution.fit(FirstPrComps,4,...
%         'Start','randSample','Replicates',5,'Option',statset('Display','final')); % 4 clusters
%     fourclus=1;
% catch
%     gmm_fit_clusters = gmdistribution.fit(FirstPrComps,3,...
%         'Start','randSample','Replicates',5,'Option',statset('Display','final')); % 3 clusters
%     fourclus=0;
% end
%
% % Look at clusters mu and sigma
% % gmm_fit_clusters.mu(1,:)
% % gmm_fit_clusters.Sigma(:,:,1)
%
% % plot cluster contours
% hold on
% clusgmmfith = ezcontour(@(x,y)pdf(gmm_fit_clusters,[x y]),[min(PrComps(:,1))-1 max(PrComps(:,1))+1],[min(PrComps(:,2))-1 max(PrComps(:,2))+1]);
% % ezsurf(@(x,y)pdf(gaussmixmodel_fit,[x y]),[min(PrComps(:,1))-1 max(PrComps(:,1))+1],[min(PrComps(:,2))-1 max(PrComps(:,2))+1]);
%
% % cluster and find posterior probability
% %[PCAclusidx,nlogl,P,logpdf,M] = cluster(gmm_fit_clusters,FirstPrComps);
% % actually use posterior value to
% ClusPr = posterior(gmm_fit_clusters,[PrComps(:,1),PrComps(:,2)]);
% % add a 5th column for junk
% ClusPr(max(ClusPr,[],2)<0.95,5)=1;
% % classify into clusters (or junk)
% [~,PCAclusidx] = max(ClusPr,[],2);
% maxClusPr = max(ClusPr(:,1:4),[],2);
% tabulate(PCAclusidx);
%
% gscatter(FirstPrComps(:,1), FirstPrComps(:,2), PCAclusidx)
% %print index number
% text(PrComps(:,1), PrComps(:,2),num2str(rot90(size(PrComps(:,1),1):-1:1)))
% %and print proba value
% text(PrComps(:,1), PrComps(:,2)-0.5,num2str(floor(maxClusPr*100)))
%
% hold off
% title('PC1 vs PC2, 4 clusters');
%
% %%
% % cluster1 = FirstPrComps(PCAclusidx == 1,:);
% % cluster2 = FirstPrComps(PCAclusidx == 2,:);
% % cluster3 = FirstPrComps(PCAclusidx == 3,:);
% % if fourclus && sum(PCAclusidx == 4)
% %     cluster4 = FirstPrComps(PCAclusidx == 4,:);
% % end
% % delete(clusgmmfith);
% % gaussfith1 = scatter(cluster1(:,1),cluster1(:,2),80,'b.');
% % gaussfith2 = scatter(cluster2(:,1),cluster2(:,2),80,'r.');
% % gaussfith3 = scatter(cluster3(:,1),cluster3(:,2),80,'g.');
% % if fourclus
% %     gaussfith4 = scatter(cluster4(:,1),cluster4(:,2),80,'c.');
% %     legend([gaussfith1 gaussfith2 gaussfith3 gaussfith4],'Cluster 1','Cluster 2','Cluster 3','Cluster 4','Location','NW')
% % else
% %     legend([gaussfith1 gaussfith2 gaussfith3],'Cluster 1','Cluster 2','Cluster 3','Location','NW')
% % end
%
% % % in 2D, with variance multiplier
% % figure
% % plot(PrComps(:,1).*latent(1),PrComps(:,2).*latent(2),'.b')
%
% %% in 3D
% figure
% scatter3(PrComps(:,1),PrComps(:,2),PrComps(:,3),40,'b.')
% xlabel('PC 1'); ylabel('PC 2'); zlabel('PC 3')
% hold on
% % outline just seeds
% scatter3(PrComps(wavedropseed,1),PrComps(wavedropseed,2),PrComps(wavedropseed,3),80,'r','filled')
% scatter3(PrComps(waveburstseed,1),PrComps(waveburstseed,2),PrComps(waveburstseed,3),80,'g','filled')
% scatter3(PrComps(waveflatseed,1),PrComps(waveflatseed,2),PrComps(waveflatseed,3),80,'k','filled')
% % outline bests
% scatter3(PrComps([6,11,21,24,29,32],1),PrComps([6,11,21,24,29,32],2),PrComps([6,11,21,24,29,32],3),80,'r','filled')
% scatter3(PrComps([7,14,17,25,27],1),PrComps([7,14,17,25,27],2),PrComps([7,14,17,25,27],3),80,'g','filled')
%
% %% plot sdfs from each cluster
% clus1=find(PCAclusidx==1); clus2=find(PCAclusidx==2);
% clus3=find(PCAclusidx==3); clus4=find(PCAclusidx==4);
% clus5=find(PCAclusidx==5);
%
% figure('name','cluster1')
% subplotdim=[ceil(sqrt(numel(clus1))),ceil(sqrt(numel(clus1)))];
% for sacplot=1:length(clus1)
%     subplot(subplotdim(1),subplotdim(2),sacplot)
%     plot(bnorm_sacresps(clus1(sacplot),:));
%     ylim=get(gca,'ylim');
%     set(gca,'ylim',[min(ylim(1),0) ylim(2)]);
%     text(10,ylim(2)-0.1,['sacplot ' num2str(clus1(sacplot))]);
% end
% figure('name','cluster2')
% subplotdim=[ceil(sqrt(numel(clus2))),ceil(sqrt(numel(clus2)))];
% for sacplot=1:length(clus2)
%     subplot(subplotdim(1),subplotdim(2),sacplot)
%     plot(bnorm_sacresps(clus2(sacplot),:));
%     ylim=get(gca,'ylim');
%     set(gca,'ylim',[min(ylim(1),0) ylim(2)]);
%     text(10,ylim(2)-0.1,['sacplot ' num2str(clus2(sacplot))]);
% end
% figure('name','cluster3')
% subplotdim=[ceil(sqrt(numel(clus3))),ceil(sqrt(numel(clus3)))];
% for sacplot=1:length(clus3)
%     subplot(subplotdim(1),subplotdim(2),sacplot)
%     plot(bnorm_sacresps(clus3(sacplot),:));
%     ylim=get(gca,'ylim');
%     set(gca,'ylim',[min(ylim(1),0) ylim(2)]);
%     text(10,ylim(2)-0.1,['sacplot ' num2str(clus3(sacplot))]);
% end
% figure('name','cluster4')
% subplotdim=[ceil(sqrt(numel(clus4))),ceil(sqrt(numel(clus4)))];
% for sacplot=1:length(clus4)
%     subplot(subplotdim(1),subplotdim(2),sacplot)
%     plot(bnorm_sacresps(clus4(sacplot),:));
%     ylim=get(gca,'ylim');
%     set(gca,'ylim',[min(ylim(1),0) ylim(2)]);
%     text(10,ylim(2)-0.1,['sacplot ' num2str(clus4(sacplot))]);
% end
% if fourclus && sum(clus5)
%     figure('name','cluster5')
%     subplotdim=[ceil(sqrt(numel(clus5))),ceil(sqrt(numel(clus5)))];
%     for sacplot=1:length(clus5)
%         subplot(subplotdim(1),subplotdim(2),sacplot)
%         plot(bnorm_sacresps(clus5(sacplot),:));
%         ylim=get(gca,'ylim');
%         set(gca,'ylim',[min(ylim(1),0) ylim(2)]);
%         text(10,ylim(2)-0.1,['sacplot ' num2str(clus5(sacplot))]);
%     end
% end

%% Hierarchical clustering
pairw_dist = pdist(rnorm_sacresps,'cityblock'); %'cityblock'
squareform(pairw_dist);
hc_links = linkage(pairw_dist,'average'); %'average'
dendrogram(hc_links);
%verify inconsistencies and dissimilarities
hc_inco = inconsistent(hc_links);
hc_dissim = cophenet(hc_links,pairw_dist);
%cluster by setting inconsistency coefficient threshold
inc_coef_th=1.15;
hc_clus = cluster(hc_links,'cutoff',inc_coef_th);
% or define number of cluster wanted
hc_clus = cluster(hc_links,'maxclust',12);

% plot HC clusters
for hclus=1:12
    figure('name',['cluster' num2str(hclus)]);
    clusn=find(hc_clus==hclus);
    subplotdim=[ceil(sqrt(numel(clusn))),ceil(sqrt(numel(clusn)))];
    for sacplot=1:length(clusn)
        subplot(subplotdim(1),subplotdim(2),sacplot)
        plot(rnorm_sacresps(clusn(sacplot),:));
        ylim=get(gca,'ylim');
        set(gca,'ylim',[min(ylim(1),0) ylim(2)]);
        text(10,ylim(2)-0.1,['sacplot ' num2str(clusn(sacplot))]);
    end
end


%% colors for population plots
%     figure(1);
cc=lines(size(allgsalignmnt,1)); % one color per file
if size(cc,1)==8
    cc(8,:)=[0 0.75 0];
end

%% prealloc compile data
compgssdf=struct('clus',{'rampfallclus','sacburstclus','rewrampclus'},...
    'align',struct('sac',struct('NSStrial',nan(1,1301),'CStrial',nan(1,1301),'NCStrial',nan(1,1301)),...
    'tgt',struct('NSStrial',nan(1,901),'CStrial',nan(1,901),'NCStrial',nan(1,901)),...
    'ssd',struct('LMCS_NSStrial',nan(1,1301),'LMNCS_NSStrial',nan(1,1301),'CStrial',nan(1,1301),'NCStrial',nan(1,1301)),...
    'corsac',struct('NSStrial',nan(1,1001),'CStrial',nan(1,1001),'NCStrial',nan(1,1001)),...
    'rew',struct('NSStrial',nan(1,1001),'CStrial',nan(1,1001),'NCStrial',nan(1,1001))));

trialtype={'NSStrial','CStrial','NCStrial'};
ssdtrialtype={'LMCS_NSStrial','LMNCS_NSStrial','CStrial','NCStrial'};

% data span
sac_startstop=[900 400];
tgt_startstop=[200 700];
ssd_startstop=[800 500];
corsac_startstop=[800 200];
rew_startstop=[800 200];

%% separate data by cluster
%% cb cx cluster
% clusgsndata{1}=allgsndata(hc_clus==4 | hc_clus==10,:);
% clusgsndata{2}=allgsndata(hc_clus==2,:);
% clusgsndata{3}=allgsndata(hc_clus==6,:);
% 
% clussblmean{1}=bslresp_mean(hc_clus==4 | hc_clus==10);
% clussblmean{2}=bslresp_mean(hc_clus==2);
% clussblmean{3}=bslresp_mean(hc_clus==6);
% 
% clussbslresp_sd{1}=bslresp_sd(hc_clus==4 | hc_clus==10);
% clussbslresp_sd{2}=bslresp_sd(hc_clus==2);
% clussbslresp_sd{3}=bslresp_sd(hc_clus==6);
% 
% clusprefdir{1}=allgsprefdir(hc_clus==4 | hc_clus==10,:);
% clusprefdir{2}=allgsprefdir(hc_clus==2,:);
% clusprefdir{3}=allgsprefdir(hc_clus==6,:);
% 
% clusssds{1}=allgsssds(hc_clus==4 | hc_clus==10);
% clusssds{2}=allgsssds(hc_clus==2);
% clusssds{3}=allgsssds(hc_clus==6);

%% cDN clusters
clusgsndata{1}=allgsndata(hc_clus==6 | hc_clus==8 | hc_clus==11,:);
clusgsndata{2}=allgsndata(hc_clus==9,:);
clusgsndata{3}=allgsndata(hc_clus==10,:);

clussblmean{1}=bslresp_mean(hc_clus==6 | hc_clus==8 | hc_clus==11);
clussblmean{2}=bslresp_mean(hc_clus==9);
clussblmean{3}=bslresp_mean(hc_clus==10);

clussbslresp_sd{1}=bslresp_sd(hc_clus==6 | hc_clus==8 | hc_clus==11);
clussbslresp_sd{2}=bslresp_sd(hc_clus==9);
clussbslresp_sd{3}=bslresp_sd(hc_clus==10);

clusprefdir{1}=allgsprefdir(hc_clus==6 | hc_clus==8 | hc_clus==11,:);
clusprefdir{2}=allgsprefdir(hc_clus==9,:);
clusprefdir{3}=allgsprefdir(hc_clus==10,:);

clusssds{1}=allgsssds(hc_clus==6 | hc_clus==8 | hc_clus==11);
clusssds{2}=allgsssds(hc_clus==9);
clusssds{3}=allgsssds(hc_clus==10);

% allgsprevssd
% allgssacdelay

% compgssdf{1}=nan(size(allgsalignmnt,1),1301,3);
% compgssdf{2}=nan(size(allgsalignmnt,1),901,3);
% compgssdf{3}=nan(size(allgsalignmnt,1),1301,4);
% compgssdf{4}=nan(size(allgsalignmnt,1),1001,3);
% compgssdf{5}=nan(size(allgsalignmnt,1),1001,3);

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
        
        if ~isempty(gsdata) && clussbslresp_sd{clusnum}(gsd)~=0
            for sacalg=1:3
                try
                    rasters=gsdata(sacalg).rast;
                    if prefdironly
                        rasters=rasters(clusprefdir{clusnum}{gsd,1}{sacalg},:); % sorted by 3 structures containing logical indices
                    end
                    alignmtt=gsdata(sacalg).alignt;
                    start=alignmtt-sac_startstop(1)-half_sixsig; stop=alignmtt+sac_startstop(2)+half_sixsig;
                    [sdf, convrasters, convrastsem]=conv_raster(rasters,conv_sigma,start,stop);
                    
                    %normalize sdf by baseline activity
                    normsdf=(sdf-clussblmean{clusnum}(gsd))./clussbslresp_sd{clusnum}(gsd);
                    
                    %% plots
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
        
        if ~isempty(gsdata) && clussbslresp_sd{clusnum}(gsd)~=0
            for tgtalg=1:3
                try
                    rasters=gsdata(tgtalg).rast;
                    if prefdironly
                        rasters=rasters(clusprefdir{clusnum}{gsd,1}{tgtalg},:); % sorted by 3 structures containing logical indices
                    end
                    alignmtt=gsdata(tgtalg).alignt;
                    start=alignmtt-tgt_startstop(1)-half_sixsig; stop=alignmtt+tgt_startstop(2)+half_sixsig;
                    [sdf, convrasters, convrastsem]=conv_raster(rasters,conv_sigma,start,stop);
                    
                    %normalize sdf by baseline activity
                    normsdf=(sdf-clussblmean{clusnum}(gsd))./clussbslresp_sd{clusnum}(gsd);
                    
                    %% plots
                    %                              figure(1)
                    %                              hold off
                    %                              patch([1:length(normsdf),fliplr(1:length(normsdf))],[normsdf-convrastsem,fliplr(normsdf+convrastsem)],'k','EdgeColor','none','FaceAlpha',0.1);
                    %                              hold on
                    %                              %plot sdf
                    %                              plot(normsdf,'Color','b','LineWidth',1.8);
                    %                              set(gca,'xtick',[1:100:1301],'xticklabel',[-800:100:500])
                    %                              close(gcf)
                    
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
                    [sdf, convrasters, convrastsem]=conv_raster(rasters,conv_sigma,start,stop);
                    
                    %normalize sdf by baseline activity
                    normsdf=(sdf-clussblmean{clusnum}(gsd))./clussbslresp_sd{clusnum}(gsd);
                    
                    %% plots
                    %          figure(1)
                    %          hold off
                    %          patch([1:length(sdf),fliplr(1:length(sdf))],[sdf-convrastsem,fliplr(sdf+convrastsem)],'k','EdgeColor','none','FaceAlpha',0.1);
                    %          hold on
                    %          %plot sdf
                    %          plot(sdf,'Color','b','LineWidth',1.8);
                    %          set(gca,'xtick',[1:100:1301],'xticklabel',[-800:100:500])
                    %          close(gcf)
                    
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
                    [sdf, convrasters, convrastsem]=conv_raster(rasters,conv_sigma,start,stop);
                    
                    %normalize sdf by baseline activity
                    normsdf=(sdf-clussblmean{clusnum}(gsd))./clussbslresp_sd{clusnum}(gsd);
                    
                    %% plots
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
                    [sdf, convrasters, convrastsem]=conv_raster(rasters,conv_sigma,start,stop);
                    
                    %normalize sdf by baseline activity
                    normsdf=(sdf-clussblmean{clusnum}(gsd))./clussbslresp_sd{clusnum}(gsd);
                    
                    %% plots
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

%% compute and plot population activity (and ci) by cluster

trialtype_sdf={'NSStrial_popsdf','CStrial_popsdf','NCStrial_popsdf'};
trialtype_ci={'NSStrial_popci','CStrial_popci','NCStrial_popci'};
ssdtrialtype_sdf={'LMCS_NSStrial_popsdf','LMNCS_NSStrial_popsdf','CStrial_popsdf','NCStrial_popsdf'};
ssdtrialtype_ci={'LMCS_NSStrial_popci','LMNCS_NSStrial_popci','CStrial_popci','NCStrial_popci'};

lineStyles = linspecer(4);
[sacfig,tgtfig,ssdfig,corsacfig,rewfig]=deal(nan(1,3)); %handles for figures

for clusnum=1:3
    
    %% saccade alignment plot
    sacfig(clusnum)=figure('name',['Cluster' num2str(clusnum) ' saccade plots']);
    hold on;
    
    for sacpop=1:3
        [popsdf, compgssdf(1,clusnum).align.sac.(trialtype_sdf{sacpop})]=deal(nanmean(compgssdf(1,clusnum).align.sac.(trialtype{sacpop})));
        [popci, compgssdf(1,clusnum).align.sac.(trialtype_ci{sacpop})]=deal(nanstd(compgssdf(1,clusnum).align.sac.(trialtype{sacpop}))/...
            sqrt(size(compgssdf(1,clusnum).align.sac.(trialtype{sacpop}),1)));
        
        patch([1:length(popci),fliplr(1:length(popci))],...
            [popsdf-popci,fliplr(popsdf+popci)],...
            lineStyles(sacpop,:),'EdgeColor','none','FaceAlpha',0.1);
        hold on;
        lineh(sacpop)=plot(popsdf,'LineWidth',2,'color',lineStyles(sacpop,:));
    end
    currylim=get(gca,'ylim');
    patch([sac_startstop(1)-2:sac_startstop(1)+2 fliplr(sac_startstop(1)-2:sac_startstop(1)+2)], ...
        reshape(repmat([currylim(1) currylim(2)],5,1),1,numel(currylim)*5), ...
        [0 0 0],'EdgeColor','none','FaceAlpha',0.5);
    hold off;
    %% beautify plot
    set(gca,'XTick',[0:100:(sac_startstop(2)+sac_startstop(1))]);
    set(gca,'XTickLabel',[-sac_startstop(1):100:sac_startstop(2)]);
    axis(gca,'tight'); box off;
    set(gca,'Color','white','TickDir','out','FontName','Cambria','FontSize',10);
    hxlabel=xlabel(gca,'Time (ms)','FontName','Cambria','FontSize',10);
    hylabel=ylabel(gca,'Firing rate (z-score)','FontName','Cambria','FontSize',10);
    legh=legend(lineh(1:3),{'No Stop Signal','Stop Signal: Cancelled', 'Stop Signal: Non Cancelled'});
    set(legh,'Interpreter','none','Location','SouthWest','Box','off','LineWidth',1.5,'FontName','Cambria','FontSize',9); % Interpreter prevents underscores turning character into subscript
    title(['Cluster' num2str(clusnum) ' Aligned to saccade'],'FontName','Cambria','FontSize',15);
    
    %% target alignment plot
    tgtfig(clusnum)=figure('name',['Cluster' num2str(clusnum) ' target plots']);
    hold on;
    
    for tgtpop=1:3
        [popsdf, compgssdf(1,clusnum).align.tgt.(trialtype_sdf{tgtpop})]=deal(nanmean(compgssdf(1,clusnum).align.tgt.(trialtype{tgtpop})));
        [popci, compgssdf(1,clusnum).align.tgt.(trialtype_ci{tgtpop})]=deal(nanstd(compgssdf(1,clusnum).align.tgt.(trialtype{tgtpop}))/...
            sqrt(size(compgssdf(1,clusnum).align.tgt.(trialtype{tgtpop}),1)));
        
        patch([1:length(popci),fliplr(1:length(popci))],...
            [popsdf-popci,fliplr(popsdf+popci)],...
            lineStyles(tgtpop,:),'EdgeColor','none','FaceAlpha',0.1);
        hold on;
        lineh(tgtpop)=plot(popsdf,'LineWidth',2,'color',lineStyles(tgtpop,:));
    end
    currylim=get(gca,'ylim');
    patch([tgt_startstop(1)-2:tgt_startstop(1)+2 fliplr(tgt_startstop(1)-2:tgt_startstop(1)+2)], ...
        reshape(repmat([currylim(1) currylim(2)],5,1),1,numel(currylim)*5), ...
        [0 0 0],'EdgeColor','none','FaceAlpha',0.5);
    hold off;
    %% beautify plot
    set(gca,'XTick',[0:100:(tgt_startstop(2)+tgt_startstop(1))]);
    set(gca,'XTickLabel',[-tgt_startstop(1):100:tgt_startstop(2)]);
    axis(gca,'tight'); box off;
    set(gca,'Color','white','TickDir','out','FontName','Cambria','FontSize',10);
    hxlabel=xlabel(gca,'Time (ms)','FontName','Cambria','FontSize',10);
    hylabel=ylabel(gca,'Firing rate (z-score)','FontName','Cambria','FontSize',10);
    legh=legend(lineh(1:3),{'No Stop Signal','Stop Signal: Cancelled', 'Stop Signal: Non Cancelled'});
    set(legh,'Interpreter','none','Location','NorthEast','Box', 'off','LineWidth',1.5,'FontName','Cambria','FontSize',9); % Interpreter prevents underscores turning character into subscript
    title(['Cluster' num2str(clusnum) ' Aligned to target'],'FontName','Cambria','FontSize',15);
    
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
        
        patch([1:length(popci),fliplr(1:length(popci))],...
            [popsdf-popci,fliplr(popsdf+popci)],...
            lineStyles(ssdpop,:),'EdgeColor','none','FaceAlpha',0.1);
        hold on;
        lineh(ssdpop)=plot(popsdf,'LineWidth',2,'color',lineStyles(ssdpop,:));
    end
    currylim=get(gca,'ylim');
    patch([ssd_startstop(1)-2:ssd_startstop(1)+2 fliplr(ssd_startstop(1)-2:ssd_startstop(1)+2)], ...
        reshape(repmat([currylim(1) currylim(2)],5,1),1,numel(currylim)*5), ...
        [0 0 0],'EdgeColor','none','FaceAlpha',0.5);
    hold off;
    %% beautify plot
    set(gca,'XTick',[0:100:(ssd_startstop(2)+ssd_startstop(1))]);
    set(gca,'XTickLabel',[-ssd_startstop(1):100:ssd_startstop(2)]);
    axis(gca,'tight'); box off;
    set(gca,'Color','white','TickDir','out','FontName','Cambria','FontSize',10);
    hxlabel=xlabel(gca,'Time (ms)','FontName','Cambria','FontSize',10);
    hylabel=ylabel(gca,'Firing rate (z-score)','FontName','Cambria','FontSize',10);
    legh=legend(lineh(1:2:3),{'Latency Matched No Stop Signal','Stop Signal: Cancelled'});
    set(legh,'Interpreter','none','Box', 'off','LineWidth',1.5,'FontName','Cambria','FontSize',9); % Interpreter prevents underscores turning character into subscript
    title(['Cluster' num2str(clusnum) ' NSS CS Aligned to ssd'],'FontName','Cambria','FontSize',15);
    
    %2nd plots
    subplot(1,2,2)
    
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
        else
            patch([1:length(popci),fliplr(1:length(popci))],...
                [popsdf-popci,fliplr(popsdf+popci)],...
                lineStyles(ssdpop,:),'EdgeColor','none','FaceAlpha',0.1);
            lineh(ssdpop)=plot(popsdf,'LineWidth',2,'color',lineStyles(ssdpop,:));
        end
    end
    currylim=get(gca,'ylim');
    patch([ssd_startstop(1)-2:ssd_startstop(1)+2 fliplr(ssd_startstop(1)-2:ssd_startstop(1)+2)], ...
        reshape(repmat([currylim(1) currylim(2)],5,1),1,numel(currylim)*5), ...
        [0 0 0],'EdgeColor','none','FaceAlpha',0.5);
    hold off;
    %% beautify plot
    set(gca,'XTick',[0:100:(ssd_startstop(2)+ssd_startstop(1))]);
    set(gca,'XTickLabel',[-ssd_startstop(1):100:ssd_startstop(2)]);
    axis(gca,'tight'); box off;
    set(gca,'Color','white','TickDir','out','FontName','Cambria','FontSize',10);
    hxlabel=xlabel(gca,'Time (ms)','FontName','Cambria','FontSize',10);
    hylabel=ylabel(gca,'Firing rate (z-score)','FontName','Cambria','FontSize',10);
    legh=legend(lineh(2:2:4),{'Latency Matched No Stop Signal','Stop Signal: Non Cancelled'});
    set(legh,'Interpreter','none','Box', 'off','LineWidth',1.5,'FontName','Cambria','FontSize',9); % Interpreter prevents underscores turning character into subscript
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
    set(legh,'Interpreter','none','Location','NorthEast','Box', 'off','LineWidth',1.5,'FontName','Cambria','FontSize',9); % Interpreter prevents underscores turning character into subscript
    title(['Cluster' num2str(clusnum) ' Aligned to corrective saccade'],'FontName','Cambria','FontSize',15);
    
    
    %% reward alignment plots
    
    rewfig(clusnum)=figure('name',['Cluster' num2str(clusnum) ' reward plots']);
    hold on;
    
    for rewpop=1:2
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
    legh=legend(lineh(1:2),{'No Stop Signal','Stop Signal: Cancelled'});
    set(legh,'Interpreter','none','Location','NorthEast','Box', 'off','LineWidth',1.5,'FontName','Cambria','FontSize',9); % Interpreter prevents underscores turning character into subscript
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











