%% Unsupervised hierarchical clustering on PCA
function [clusterIdx,clusterSort,zTraces]=UHC_PCA(traces,SNRs,numClus)

%% zscore traces
zTraces=zscore(traces,[],2);

%% Or use PCA
% https://www.sciencedirect.com/science/article/pii/S0896627317303434
[coeffs,~,~,~,explained] = pca(traces(SNRs>5,:)'); %or zTraces  

% ninetyfiveVariancePoint=find(cumsum(explained)>=95,1)-1;
%Euclidean distance
eDistance = pdist(coeffs(:,1:4)); %'cityblock'

% hierarchical clustering tree
clustTree= linkage(eDistance,'complete'); %'average'
% leafOrder = optimalleaforder(clustTree,eDistance); 
%,'Criteria','group','Transformation','inverse');
% clusters
if ~exist('numClus','var')
    numClus=3;
end
clusterIdx = cluster(clustTree,'maxclust',numClus); %'criterion','distance','cutoff',2.5*10^4

[proportions,clusterIDs]=hist(clusterIdx,unique(clusterIdx))

%% now classify all traces into those clusters
[clusterSortedIdx,clusterSort]=sort(clusterIdx);

% sort clusters as 3/4/1/2
clusterSort=[clusterSort(clusterSortedIdx==3);clusterSort(clusterSortedIdx==1);...
    clusterSort(clusterSortedIdx==4);clusterSort(clusterSortedIdx==2)];

distanceMatrix=pdist(zTraces(clusterSort,:),'correlation');

% clusterIDs=[1;4;2;3]; %[2;1;3]; %[3;4;1;2]

%% summary figure
figure;
%order by cluster to plot heatmap
subplot(3,numClus*2,[1:numClus,numClus*2+1:numClus*3]);
% imagesc(squareform(distanceMatrix)); %colormap('jet') % for redbluecmap: install  Bioinformatics Toolbox
heatmap(squareform(distanceMatrix),'ColorMethod','median');
colorbar('northoutside')
if size(zTraces,2) == 500 %short
    tickVals=100:100:500;
    tickLabelVals=[-200 -100 0 100 200];
    patchLoc=[299 300 300 299];
elseif size(zTraces,2) == 700 %long
    tickVals=100:100:700;
    tickLabelVals=[-400 -300 -200 -100 0 100 200];
    patchLoc=[499 500 500 499];
end
cmap=lines; cmap=[0,0,0;cmap];
clusterProfile=nan(numClus,size(zTraces,2));
for clusNum = 1:numClus
    subplot(3,numClus*2,numClus*4+clusNum);
    clust = clusterIdx==clusterIDs(clusNum);
    if size(zTraces(clust,:),1)>1
        clusterProfile(clusNum,:)=nanmean(zTraces(clust,:));
    else
        clusterProfile(clusNum,:)=zTraces(clust,:);
    end
    plot(clusterProfile(clusNum,:),'Color',cmap(clusNum+1,:),'Linewidth',2);
    set(gca,'xtick', tickVals, 'xticklabel',tickLabelVals)
    yLims=get(gca,'YLim');
    %     set(gca,'YLim',[0 max(yLims)]);
    patch(patchLoc,[min(yLims) max(yLims) max(yLims) min(yLims)],...
        cmap(clusNum,:),'FaceAlpha',0.5);
    text(1,round(max(yLims)-1),['n=' num2str(size(zTraces(clust,:),1))]);
    axis tight
end

subplot(2,numClus*2,[numClus+1:numClus*2,numClus*3+1:numClus*4]);
imagesc(zTraces(clusterSort,:)); axis tight; yLims=get(gca,'YLim');
patch(patchLoc,[min(yLims) max(yLims) max(yLims) min(yLims)],...
    cmap(clusNum,:),'FaceAlpha',0.5);
    set(gca,'xtick', tickVals, 'xticklabel',tickLabelVals)
colorbar;
