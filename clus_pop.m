function [clusidx,clustypes,clusavwf]=clus_pop(sacresps,bnorm_sacresps,rnorm_sacresps,method,option)

if nargin<4
    method='hclus';
    option='round1';
end
if nargin<5
    option='round1';
end

%% smooth sdfs
sm_bnorm_sacresps=bnorm_sacresps;
sm_rnorm_sacresps=rnorm_sacresps;
for resp=1:size(bnorm_sacresps,1)
    sm_bnorm_sacresps(resp,:) = gauss_filtconv(bnorm_sacresps(resp,:),50);
    sm_rnorm_sacresps(resp,:) = gauss_filtconv(rnorm_sacresps(resp,:),50);
end

%evaluate optimal number of clusters
davies_index = zeros(1,15);
davies_index(1) = 999;
for num_clusters = 2:15
    [cluster_ids,cluster_centroids,cluster_distances] = kmeans(sm_rnorm_sacresps,num_clusters,'Distance','sqEuclidean');
    
    pairwise_cluster_distances = squareform(pdist(cluster_centroids,'Euclidean'));
    cluster_populations = arrayfun(@(x) numel(find(cluster_ids == x)),1:num_clusters);
    intracluster_distances = cluster_distances./cluster_populations';
    
    davies_per_cluster = zeros(1,num_clusters);
    for clusnb = 1:num_clusters
        davies_temp = zeros(1,num_clusters);
        for j= setxor(1:num_clusters,clusnb)
            davies_temp(j) = (intracluster_distances(clusnb)+intracluster_distances(j))./pairwise_cluster_distances(j,clusnb);
        end
        davies_per_cluster(clusnb) = max(davies_temp);
    end
    davies_index(num_clusters) = 1./num_clusters*sum(davies_per_cluster);
end
% figure; plot(davies_index(2:end),'r','linewidth',2)
opti_clusnm=find(davies_index(2:end)==max(davies_index(2:end)))+1;

%% activity profiles (used for seeds)
midrange=size(bnorm_sacresps,2)/2;

% Make seed that represent midrange drop / midrange burst / ramp to end
midrangeseeds=cellfun(@(x) mean(x(1,midrange-150:midrange-50))-mean(x(1,midrange+50:midrange+150)), mat2cell(bnorm_sacresps,ones(size(bnorm_sacresps,1),1)));
peakseeds=cellfun(@(x) (mean(x(1,midrange+50:midrange+100))-mean(x(1,midrange-150:midrange-50)))+...
    (mean(x(1,midrange+50:midrange+100))-mean(x(1,midrange+100:midrange+200))), mat2cell(bnorm_sacresps,ones(size(bnorm_sacresps,1),1)));
outerrangeseeds=cellfun(@(x) mean(x(1,length(x)-150:length(x)-1))-mean(x(1,1:150)), mat2cell(bnorm_sacresps,ones(size(bnorm_sacresps,1),1)));

% keep 10 highest seed values
[~,mrseeds_vals_idx]=sort(midrangeseeds);
[~,pkseeds_vals_idx]=sort(peakseeds);
[~,orseeds_vals_idx]=sort(outerrangeseeds);
drop_seeds_vals=bnorm_sacresps(mrseeds_vals_idx(end-10:end),:);
burst_seeds_vals=bnorm_sacresps(pkseeds_vals_idx(end-10:end),:);
rampatw_seeds_vals=bnorm_sacresps(orseeds_vals_idx(end-10:end),:);

seeds=[nanmean(drop_seeds_vals);...
    nanmean(burst_seeds_vals);...
    nanmean(rampatw_seeds_vals)];

% create parameters for fit template function
xfit_vals=linspace(1,size(seeds,2),size(seeds,2));
drop_seed_polyf = polyfit(xfit_vals,seeds(1,:),5);
% drop_yfit_vals=polyval(drop_seed_polyf,xfit_vals);
burst_seed_polyf = polyfit(xfit_vals,seeds(2,:),5);
% burst_yfit_vals=polyval(burst_seed_polyf,xfit_vals);
rampatw_seed_polyf = polyfit(xfit_vals,seeds(3,:),5);
% rampatw_yfit_vals=polyval(rampatw_seed_polyf,xfit_vals);

% for i=1:20
%     figure; plot(bnorm_sacresps(i,:));
%     text(10,min(get(gca,'ylim'))+1,['var ' num2str(var(bnorm_sacresps(i,:)))]);
% end

%% plot seeds
% figure
% plot(seeds(1,:))
% hold on
% plot(seeds(2,:),'r')
% plot(seeds(3,:),'g')


%%%%%%%%%%%%%%%%%%%%%%
%% clusterization
%%%%%%%%%%%%%%%%%%%%%%

if strcmp(method,'kmean')
    %% k-means clustering
    
    %% calculate k-means
    % at random
    % [IDX,C,sumd,D]=kmeans(sacresps,3,'dist','city','display','iter');
    % seeded
    [kMidx,kMeansClus,sumd,D]=kmeans(bnorm_sacresps,opti_clusnm,'dist','city','start',seeds,'display','iter');
    %try [cluster_ids] = kmeans(bnorm_sacresps,opti_clusnm,'Distance','sqEuclidean');
    %% plot means
    figure
    plot(kMeansClus(1,:),'b')
    hold on
    plot(kMeansClus(2,:),'r')
    plot(kMeansClus(3,:),'g')
    
    %% k-means scatterplot
    figure
    scatter3(D(kMidx==1,1),D(kMidx==1,2),D(kMidx==1,3),40,'b.')
    hold on
    scatter3(D(kMidx==2,1),D(kMidx==2,2),D(kMidx==2,3),40,'r.')
    scatter3(D(kMidx==3,1),D(kMidx==3,2),D(kMidx==3,3),40,'g.')
    
    % Fit Gaussian Mixture Model using the k-means centers as the initial conditions
    % We only have mean initial conditions from the k-means algorithm, so we
    % can specify some arbitrary initial variance and mixture weights.
    % gmInitialVariance = 0.1;
    % initialSigma = cat(3,gmInitialVariance,gmInitialVariance,gmInitialVariance);
    % % Initial weights are set at 50%
    % initialWeights = [0.5 0.5 0.5];
    % % Initial condition structure for the gmdistribution.fit function
    % S.mu = kMeansClus;
    % S.Sigma = initialSigma;
    % S.PComponents = initialWeights;
    
    %% plot sdfs from each k-means cluster
    clus1=find(kMidx==1); clus2=find(kMidx==2); clus3=find(kMidx==3);
    figure('name','cluster1')
    subplotdim=[ceil(length(clus1)/2)-(2*floor(length(clus1)/10)),2+floor(length(clus1)/10)];
    for sacplot=1:length(clus1)
        subplot(subplotdim(1),subplotdim(2),sacplot)
        plot(bnorm_sacresps(clus1(sacplot),:));
        ylim=get(gca,'ylim');
        set(gca,'ylim',[0 ylim(2)]);
        text(20,10,['sacplot ' num2str(clus1(sacplot))]);
    end
    figure('name','cluster2')
    subplotdim=[ceil(length(clus2)/2)-(2*floor(length(clus2)/10)),2+floor(length(clus2)/10)];
    for sacplot=1:length(clus2)
        subplot(subplotdim(1),subplotdim(2),sacplot)
        plot(bnorm_sacresps(clus2(sacplot),:));
        ylim=get(gca,'ylim');
        set(gca,'ylim',[0 ylim(2)]);
        text(20,10,['sacplot ' num2str(clus2(sacplot))]);
    end
    figure('name','cluster3')
    subplotdim=[ceil(length(clus3)/2)-(2*floor(length(clus3)/10)),2+floor(length(clus3)/10)];
    for sacplot=1:length(clus3)
        subplot(subplotdim(1),subplotdim(2),sacplot)
        plot(bnorm_sacresps(clus3(sacplot),:));
        ylim=get(gca,'ylim');
        set(gca,'ylim',[0 ylim(2)]);
        text(20,10,['sacplot ' num2str(clus3(sacplot))]);
    end
    
elseif strcmp(method,'PCA')
    
    %% PCA clustering
    [coeffs,PrComps,latent] = pca(bnorm_sacresps);
    
    % D=coeffs(:,1:8)';
    % figure
    % plot(coeffs(:,1),'.','MarkerSize',.5)
    % hold on
    % plot(coeffs(:,2),'.','MarkerSize',.5)
    % plot(coeffs(:,3),'.','MarkerSize',.5)
    % plot(coeffs(:,4),'.','MarkerSize',.5)
    
    %% scatterplot
    % in 2D
    figure
    scatter(PrComps(:,1), PrComps(:,2), 'k.');
    xlabel('PC 1'); ylabel('PC 2')
    % hold on
    % scatter(PrComps([6,11,21,24,29,32],1),PrComps([6,11,21,24,29,32],2),80,'r','filled')
    % scatter(PrComps([7,14,17,25,27],1),PrComps([7,14,17,25,27],2),80,'g','filled')
    
    %% isolate clusters
    FirstPrComps=[PrComps(:,1),PrComps(:,2)];
    try
        gmm_fit_clusters = gmdistribution.fit(FirstPrComps,4,...
            'Start','randSample','Replicates',5,'Option',statset('Display','final')); % 4 clusters
        fourclus=1;
    catch
        gmm_fit_clusters = gmdistribution.fit(FirstPrComps,3,...
            'Start','randSample','Replicates',5,'Option',statset('Display','final')); % 3 clusters
        fourclus=0;
    end
    
    % Look at clusters mu and sigma
    % gmm_fit_clusters.mu(1,:)
    % gmm_fit_clusters.Sigma(:,:,1)
    
    % plot cluster contours
    hold on
    clusgmmfith = ezcontour(@(x,y)pdf(gmm_fit_clusters,[x y]),[min(PrComps(:,1))-1 max(PrComps(:,1))+1],[min(PrComps(:,2))-1 max(PrComps(:,2))+1]);
    % ezsurf(@(x,y)pdf(gaussmixmodel_fit,[x y]),[min(PrComps(:,1))-1 max(PrComps(:,1))+1],[min(PrComps(:,2))-1 max(PrComps(:,2))+1]);
    
    % cluster and find posterior probability
    %[PCAclusidx,nlogl,P,logpdf,M] = cluster(gmm_fit_clusters,FirstPrComps);
    % actually use posterior value to
    ClusPr = posterior(gmm_fit_clusters,[PrComps(:,1),PrComps(:,2)]);
    % add a 5th column for junk
    ClusPr(max(ClusPr,[],2)<0.95,5)=1;
    % classify into clusters (or junk)
    [~,PCAclusidx] = max(ClusPr,[],2);
    maxClusPr = max(ClusPr(:,1:4),[],2);
    tabulate(PCAclusidx);
    
    gscatter(FirstPrComps(:,1), FirstPrComps(:,2), PCAclusidx)
    %print index number
    text(PrComps(:,1), PrComps(:,2),num2str(rot90(size(PrComps(:,1),1):-1:1)))
    %and print proba value
    text(PrComps(:,1), PrComps(:,2)-0.5,num2str(floor(maxClusPr*100)))
    
    hold off
    title('PC1 vs PC2, 4 clusters');
    
    %%
    % cluster1 = FirstPrComps(PCAclusidx == 1,:);
    % cluster2 = FirstPrComps(PCAclusidx == 2,:);
    % cluster3 = FirstPrComps(PCAclusidx == 3,:);
    % if fourclus && sum(PCAclusidx == 4)
    %     cluster4 = FirstPrComps(PCAclusidx == 4,:);
    % end
    % delete(clusgmmfith);
    % gaussfith1 = scatter(cluster1(:,1),cluster1(:,2),80,'b.');
    % gaussfith2 = scatter(cluster2(:,1),cluster2(:,2),80,'r.');
    % gaussfith3 = scatter(cluster3(:,1),cluster3(:,2),80,'g.');
    % if fourclus
    %     gaussfith4 = scatter(cluster4(:,1),cluster4(:,2),80,'c.');
    %     legend([gaussfith1 gaussfith2 gaussfith3 gaussfith4],'Cluster 1','Cluster 2','Cluster 3','Cluster 4','Location','NW')
    % else
    %     legend([gaussfith1 gaussfith2 gaussfith3],'Cluster 1','Cluster 2','Cluster 3','Location','NW')
    % end
    
    % % in 2D, with variance multiplier
    % figure
    % plot(PrComps(:,1).*latent(1),PrComps(:,2).*latent(2),'.b')
    
    %% in 3D
    figure
    scatter3(PrComps(:,1),PrComps(:,2),PrComps(:,3),40,'b.')
    xlabel('PC 1'); ylabel('PC 2'); zlabel('PC 3')
    hold on
    % outline just seeds
    scatter3(PrComps(wavedropseed,1),PrComps(wavedropseed,2),PrComps(wavedropseed,3),80,'r','filled')
    scatter3(PrComps(waveburstseed,1),PrComps(waveburstseed,2),PrComps(waveburstseed,3),80,'g','filled')
    scatter3(PrComps(waveflatseed,1),PrComps(waveflatseed,2),PrComps(waveflatseed,3),80,'k','filled')
    % outline bests
    scatter3(PrComps([6,11,21,24,29,32],1),PrComps([6,11,21,24,29,32],2),PrComps([6,11,21,24,29,32],3),80,'r','filled')
    scatter3(PrComps([7,14,17,25,27],1),PrComps([7,14,17,25,27],2),PrComps([7,14,17,25,27],3),80,'g','filled')
    
    %% plot sdfs from each cluster
    clus1=find(PCAclusidx==1); clus2=find(PCAclusidx==2);
    clus3=find(PCAclusidx==3); clus4=find(PCAclusidx==4);
    clus5=find(PCAclusidx==5);
    
    figure('name','cluster1')
    subplotdim=[ceil(sqrt(numel(clus1))),ceil(sqrt(numel(clus1)))];
    for sacplot=1:length(clus1)
        subplot(subplotdim(1),subplotdim(2),sacplot)
        plot(bnorm_sacresps(clus1(sacplot),:));
        ylim=get(gca,'ylim');
        set(gca,'ylim',[min(ylim(1),0) ylim(2)]);
        text(10,ylim(2)-0.1,['sacplot ' num2str(clus1(sacplot))]);
    end
    figure('name','cluster2')
    subplotdim=[ceil(sqrt(numel(clus2))),ceil(sqrt(numel(clus2)))];
    for sacplot=1:length(clus2)
        subplot(subplotdim(1),subplotdim(2),sacplot)
        plot(bnorm_sacresps(clus2(sacplot),:));
        ylim=get(gca,'ylim');
        set(gca,'ylim',[min(ylim(1),0) ylim(2)]);
        text(10,ylim(2)-0.1,['sacplot ' num2str(clus2(sacplot))]);
    end
    figure('name','cluster3')
    subplotdim=[ceil(sqrt(numel(clus3))),ceil(sqrt(numel(clus3)))];
    for sacplot=1:length(clus3)
        subplot(subplotdim(1),subplotdim(2),sacplot)
        plot(bnorm_sacresps(clus3(sacplot),:));
        ylim=get(gca,'ylim');
        set(gca,'ylim',[min(ylim(1),0) ylim(2)]);
        text(10,ylim(2)-0.1,['sacplot ' num2str(clus3(sacplot))]);
    end
    figure('name','cluster4')
    subplotdim=[ceil(sqrt(numel(clus4))),ceil(sqrt(numel(clus4)))];
    for sacplot=1:length(clus4)
        subplot(subplotdim(1),subplotdim(2),sacplot)
        plot(bnorm_sacresps(clus4(sacplot),:));
        ylim=get(gca,'ylim');
        set(gca,'ylim',[min(ylim(1),0) ylim(2)]);
        text(10,ylim(2)-0.1,['sacplot ' num2str(clus4(sacplot))]);
    end
    if fourclus && sum(clus5)
        figure('name','cluster5')
        subplotdim=[ceil(sqrt(numel(clus5))),ceil(sqrt(numel(clus5)))];
        for sacplot=1:length(clus5)
            subplot(subplotdim(1),subplotdim(2),sacplot)
            plot(bnorm_sacresps(clus5(sacplot),:));
            ylim=get(gca,'ylim');
            set(gca,'ylim',[min(ylim(1),0) ylim(2)]);
            text(10,ylim(2)-0.1,['sacplot ' num2str(clus5(sacplot))]);
        end
    end
    
elseif strcmp(method,'hclus')
    %% Hierarchical clustering
    
    pairw_dist = pdist(sm_rnorm_sacresps,'cityblock'); %'cityblock'
    squareform(pairw_dist);
    hc_links = linkage(pairw_dist,'average'); %'average'
    % dendrogram(hc_links);
    
    %verify inconsistencies and dissimilarities
    hc_inco = inconsistent(hc_links);
    hc_dissim = cophenet(hc_links,pairw_dist);
    
    %cluster by setting inconsistency coefficient threshold
    % inc_coef_th=1.15;
    % hc_clus = cluster(hc_links,'cutoff',inc_coef_th);
    % or define number of cluster wanted
    hc_clus = cluster(hc_links,'maxclust',6); %opti_clusnm
    
    % plot HC clusters
%     for hclus=1:6
%         figure('name',['cluster' num2str(hclus)]);
%         clusn=find(hc_clus==hclus);
%         subplotdim=[ceil(sqrt(numel(clusn))),ceil(sqrt(numel(clusn)))];
%         for sacplot=1:length(clusn)
%             subplot(subplotdim(1),subplotdim(2),sacplot)
%             plot(rnorm_sacresps(clusn(sacplot),:));
%             ylim=get(gca,'ylim');
%             set(gca,'ylim',[min(ylim(1),0) ylim(2)]);
%             text(10,ylim(2)-0.1,['sacplot ' num2str(clusn(sacplot))]);
%         end
%         % diffs
%         %     figure('name',['cluster' num2str(hclus) 'diffs']);
%         %     clusn=find(hc_clus==hclus);
%         %     subplotdim=[ceil(sqrt(numel(clusn))),ceil(sqrt(numel(clusn)))];
%         %     for sacplot=1:length(clusn)
%         %         subplot(subplotdim(1),subplotdim(2),sacplot)
%         %         smoothresp = gauss_filtconv(rnorm_sacresps(clusn(sacplot),:),50);
%         %         plot(diff(smoothresp));
%         %         ylim=get(gca,'ylim');
%         %         set(gca,'ylim',[min(ylim(1),0) ylim(2)]);
%         %         text(10,ylim(2)-0.1,['sacplot ' num2str(clusn(sacplot))]);
%         %     end
%     end
end

switch method
    case 'kmeans'
        clusidx=kMidx;
    case 'PCA'
        clusidx=PCAclusidx;
    case 'hclus'
        clusidx=hc_clus;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% cluster categories
%%%%%%%%%%%%%%%%%%%%%%%%%%

%% refine clusters and save average waveform

for clus=1:max(clusidx)
    %get cluster responses
    subclusidx=find(clusidx==clus);
    clusresps=sm_bnorm_sacresps(subclusidx,:);
    
    %add seeds to sub-cluster
    % clusresps=[drop_seeds_vals;burst_seeds_vals;rampatw_seeds_vals;clusresps];
    
    %     drop_seeds_vals=gauss_filtconv(mean(bnorm_sacresps(mrseeds_vals_idx(end-10:end),:)), 50);
    %     burst_seeds_vals=gauss_filtconv(mean(bnorm_sacresps(pkseeds_vals_idx(end-10:end),:)), 50);
    %     rampatw_seeds_vals=gauss_filtconv(mean(bnorm_sacresps(orseeds_vals_idx(end-10:end),:)), 50);
    %
    % find best template by minimizing shift
    varminshift=nan(size(clusresps,1),4);
    varminshift(:,1)=var(clusresps,0,2);
    for respsnm=1:size(clusresps,1)
        [~,varminshift(respsnm,2)] = fminsearch(@(shift) template_curve_match(shift,xfit_vals,clusresps(respsnm,:),drop_seed_polyf), 250);
        [~,varminshift(respsnm,3)] = fminsearch(@(shift) template_curve_match(shift,xfit_vals,clusresps(respsnm,:),burst_seed_polyf), 250);
        [~,varminshift(respsnm,4)] = fminsearch(@(shift) template_curve_match(shift,xfit_vals,clusresps(respsnm,:),rampatw_seed_polyf), 250);
    end
    
    
    if strcmp(option,'round1')
        %tag responses with variance <1 to trash/recluster pool
        clusidx(subclusidx(varminshift(:,1)<1))=-1;
        %tag other responses according to best template matching
        [~,besttempl]=min(varminshift(:,2:4),[],2);
        clusidx(subclusidx(varminshift(:,1)>1))=besttempl(varminshift(:,1)>1)+100;
    elseif strcmp(option,'round2')
        %find way to reject junk and tag others
        [~,besttempl]=min(varminshift(:,2:4),[],2);
        clusidx(subclusidx(varminshift(:,1)>1))=besttempl(varminshift(:,1)>1)+100;
    end
    

    
    %     cleanedupclus= sm_bnorm_sacresps(subclusidx([drop_minshift{:,2}]<10),:);
    %     for i=1:33
    %     figure; plot(cleanedupclus(i,:));
    %     text(10,min(get(gca,'ylim'))+1,['var ' num2str(var(cleanedupclus(i,:)))]);
    %     end
    
end

    %% define cluster type
    clustypes={'rampup','sacburst','ramp_to_reward','junk'};
    allclustypes=cell(size(clusidx));
    allclustypes(clusidx==101)=clustypes(1);
    allclustypes(clusidx==102)=clustypes(2);
    allclustypes(clusidx==103)=clustypes(3);
    allclustypes(clusidx==-1)=clustypes(4);

%% recluster discarded responses
if strcmp(option,'round1')
    [clusidx,allclustypes,clusavwf]=clus_pop(sacresps(clusidx==-1,:),bnorm_sacresps(clusidx==-1,:),rnorm_sacresps(clusidx==-1,:),method,'round2')
else
    clusavwf=[];
    return
end

%% graveyard of past subclustering attempts
% %evaluate optimal number of clusters
% davies_index = zeros(1,15);
% davies_index(1) = 999;
% for num_clusters = 2:15
%     [cluster_ids,cluster_centroids,cluster_distances] = kmeans(sm_bnorm_sacresps(trashclus,:),num_clusters,'Distance','sqEuclidean');
%
%     pairwise_cluster_distances = squareform(pdist(cluster_centroids,'Euclidean'));
%     cluster_populations = arrayfun(@(x) numel(find(cluster_ids == x)),1:num_clusters);
%     intracluster_distances = cluster_distances./cluster_populations';
%
%     davies_per_cluster = zeros(1,num_clusters);
%     for clusnb = 1:num_clusters
%         davies_temp = zeros(1,num_clusters);
%         for j= setxor(1:num_clusters,clusnb)
%             davies_temp(j) = (intracluster_distances(clusnb)+intracluster_distances(j))./pairwise_cluster_distances(j,clusnb);
%         end
%         davies_per_cluster(clusnb) = max(davies_temp);
%     end
%     davies_index(num_clusters) = 1./num_clusters*sum(davies_per_cluster);
% end
% % figure; plot(davies_index(2:end),'r','linewidth',2)
% opti_clusnm=find(davies_index(2:end)==max(davies_index(2:end)))+1;
%
% % isolate clusters using PCA
% [~,PrComps] = pca(sm_bnorm_sacresps(trashclus,:));
% FirstPrComps=[PrComps(:,1),PrComps(:,2)];
% gmm_fit_clusters = gmdistribution.fit(FirstPrComps,4,...
%     'Start','randSample','Replicates',5);
% if strcmp(lastwarn,...
%         ['Failed to converge in 100 iterations duringreplicate 4 for gmdistribution with ' num2str(4) ' components'])
%     %try again
% gmm_fit_clusters = gmdistribution.fit(FirstPrComps,4,...
%     'Start','randSample','Replicates',5);
% end
%
% % cluster and find posterior probability
% ClusPr = posterior(gmm_fit_clusters,[PrComps(:,1),PrComps(:,2)]);
% % add a n+1th column for junk
% ClusPr(max(ClusPr,[],2)<0.95,4+1)=1;
% % classify into clusters (or junk)
% [~,PCAclusidx] = max(ClusPr,[],2);
%
%
% figure
% scatter(PrComps(:,1), PrComps(:,2), 'k.');
% xlabel('PC 1'); ylabel('PC 2')
% hold on
% clusgmmfith = ezcontour(@(x,y)pdf(gmm_fit_clusters,[x y]),[min(PrComps(:,1))-1 max(PrComps(:,1))+1],[min(PrComps(:,2))-1 max(PrComps(:,2))+1]);
% %
% gscatter(FirstPrComps(:,1), FirstPrComps(:,2), PCAclusidx)
% text(PrComps(:,1), PrComps(:,2),num2str(rot90(size(PrComps(:,1),1):-1:1)))
%
%
% trashresps=sm_bnorm_sacresps(trashclus,:);
% figure;
%  plot(mean(trashresps(PCAclusidx==1,:)))
%  hold on;
%  plot(mean(trashresps(PCAclusidx==2,:)),'r')
%  plot(mean(trashresps(PCAclusidx==3,:)),'k')
%
% realjunk=trashresps(PCAclusidx==3,:);
%
%         for i=1:56
%         figure; plot(realjunk(i,:));
%         text(10,min(get(gca,'ylim'))+1,['var ' num2str(realjunk(i,:))]);
%         end
%
%
% D = mahal(gmm_fit_clusters,FirstPrComps);
% delete(clusgmmfith)
% scatter(FirstPrComps(:,1),FirstPrComps(:,2),50,D(:,1),'.')
% hb = colorbar;
% ylabel(hb,'Mahalanobis Distance to Component 1')



% % keep means
% subclusmeans=[nanmean(clusresps(PCAclusidx==1,:));...
%     nanmean(clusresps(PCAclusidx==2,:));...
%     nanmean(clusresps(PCAclusidx==3,:));...
%     nanmean(clusresps(PCAclusidx==4,:))];
%
% %% compare to reference (seeds), using mahal distance
% % needs to reduce vectors to simpler values: take max diffs, and cumsu
% subclusmeans_vals=[round(max(diff(subclusmeans,1,2),[],2).*10000) round(max(cumsum(abs(subclusmeans),2),[],2))];
%
% drop_seeds_vals=[round(max(diff(drop_seeds_vals,1,2),[],2).*10000) round(max(cumsum(abs(drop_seeds_vals),2),[],2))];
% burst_seeds_vals=[round(max(diff(burst_seeds_vals,1,2),[],2).*10000) round(max(cumsum(abs(burst_seeds_vals),2),[],2))];
% rampatw_seeds_vals=[round(max(diff(rampatw_seeds_vals,1,2),[],2).*10000) round(max(cumsum(abs(rampatw_seeds_vals),2),[],2))];
% % get mahalanobis value for each seed comparison
% mahald = [mahal(subclusmeans_vals,drop_seeds_vals) mahal(subclusmeans_vals,burst_seeds_vals) mahal(subclusmeans_vals,rampatw_seeds_vals)];
% % keep subclusters with mean mahal distance smaller than 3
% % classify the rest or send to a common pool
% resubclus=find([mean(round(mahald),2)>3 | isnan(mean(mahald,2))]);
% for reclus=1:size(resubclus,1)
%     clusidx(subclusidx(PCAclusidx==resubclus(reclus)))=-1;
% end


% find template "fitting" clusters

% Not sure about using Mahal here.    tf_clus=find(sum(mahald==min(mahald(:,1)) |  mahald==min(mahald(:,2)) | mahald==min(mahald(:,3)),2));
% Include template in clustering, get distance and residuals, classify accordingly

%% save average waveform
clusavwf=nan(length(unique(clusidx)),size(rnorm_sacresps,2));


clusavwf(clus,:)=mean(rnorm_sacresps(clusidx==clus,:));

end


