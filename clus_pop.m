%% functional clustering
function [clusterIdx,clusterSort,zTraces]=clus_pop(traces,SNRs,numClus)

if nargin==1
    SNRs=[];
end

%% zscore traces
zTraces=zscore(traces,[],2);

%% "seed" traces
% if SNRs are provided, make first round of clustering based on top 50 cells
if ~isempty(SNRs)
    [~,SeedCells_Idx]=sort(SNRs,'descend'); SeedCells_Idx=...
        SeedCells_Idx(1:round(length(SeedCells_Idx)/100*30));
    % lousyOnes=~ismember(1:size(traces,1),SeedCells_Idx);
else
    %     lousyOnes=[];
    SeedCells_Idx=1:size(traces,1);
end

seedTraces=zTraces(SeedCells_Idx,:);

%Euclidean distance
eDistance = pdist(seedTraces); %'cityblock'
% hierarchical clustering tree
clustTree= linkage(eDistance,'complete'); %'average'
% leafOrder = optimalleaforder(clustTree,eDistance); 
%,'Criteria','group','Transformation','inverse');
% clusters
if ~exist('numClus','var')
    numClus=3;
end
clusIdx = cluster(clustTree,'maxclust',numClus); %'criterion','distance','cutoff',2.5*10^4

[proportions,clusterIDs]=hist(clusIdx,unique(clusIdx))

%% find profiles
clusterProfile=nan(numClus,size(seedTraces,2));
for clusNum = 1:numClus
    clust = clusIdx==clusterIDs(clusNum);
    if size(seedTraces(clust,:),1)>1
        clusterProfile(clusNum,:)=nanmean(seedTraces(clust,:));
    else
        clusterProfile(clusNum,:)=seedTraces(clust,:);
    end
end

%% seed cluster figure
% figure;
% %order by cluster to plot heatmap
% subplot(3,numClus*2,[1:numClus,numClus*2+1:numClus*3])
% [~,clusterSort]=sort(clusIdx);
% eDistance = pdist(seedTraces(clusterSort,:)); %'cityblock'
% imagesc(squareform(eDistance)), colormap(hot); colorbar
% % figure; heatmap(squareform(eDistance));
% 
% cmap=lines; cmap=[0,0,0;cmap];
% clusterProfile=nan(numClus,size(seedTraces,2));
% for clusNum = 1:numClus
%     subplot(3,numClus*2,numClus*4+clusNum);
%     clust = clusIdx==clusterIDs(clusNum);
%     if size(seedTraces(clust,:),1)>1
%         clusterProfile(clusNum,:)=nanmean(seedTraces(clust,:));
%     else
%         clusterProfile(clusNum,:)=seedTraces(clust,:);
%     end
%     plot(clusterProfile(clusNum,:));
%     yLims=get(gca,'YLim');
%     %     set(gca,'YLim',[0 max(yLims)]);
%     patch([299 300 300 299],[min(yLims) max(yLims) max(yLims) min(yLims)],...
%         cmap(clusNum,:),'FaceAlpha',0.5);
%     text(1,round(max(yLims)-1),['n=' num2str(size(zTraces(clust,:),1))]);
% end
% 
% subplot(2,numClus*2,[numClus+1:numClus*2,numClus*3+1:numClus*4])
% imagesc(seedTraces(clusterSort,:)); colorbar;

%% now classify all traces into those clusters
profileCorrelation=corr(zTraces',clusterProfile');
[~,clusterIdx]=min(1-profileCorrelation,[],2);
[clusterSortedIdx,clusterSort]=sort(clusterIdx);
% [~,initialRank]=sort(clusterSort);
% figure; hold on
% plot(nanmean(zTraces(clusterIdx==1,:)));
% plot(nanmean(zTraces(clusterSortedIdx(initialRank)==1,:)))

% sort each cluster by trough time
% for clusNum=1:numClus
%     clusTraces=zTraces(clusterIdx==clusNum,:);
%     [~,troughSortIdx]=min(clusTraces,[],2);
%     [~,clusOrder]=sort(troughSortIdx);
%     thatClusterIdx=clusterSort(clusterSortedIdx==clusNum);
%     clusterSort(clusterSortedIdx==clusNum)=thatClusterIdx(clusOrder);
% end

% sort clusters as 3/4/1/2
clusterSort=[clusterSort(clusterSortedIdx==3);clusterSort(clusterSortedIdx==1);...
    clusterSort(clusterSortedIdx==4);clusterSort(clusterSortedIdx==2)];
% clusterIdx=[clusterSortedIdx(clusterSortedIdx==4);clusterSortedIdx(clusterSortedIdx==3);...
%     clusterSortedIdx(clusterSortedIdx==2);clusterSortedIdx(clusterSortedIdx==1)];
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
