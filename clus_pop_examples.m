%% functional clustering

function clusters=clus_pop_examples(traces,baselines)

cmap=lines; cmap=[0,0,0;cmap];

% Deneux et al 2016
% https://www.nature.com/articles/ncomms12682/figures/4
% We performed a hierarchical clustering of significant single neuron responses,
% using the similarity of temporal response profiles across neurons as a metric (Methods).
% Response traces for all sounds were smoothed using a Gaussian filter (?=31?ms).
% Before clustering, we selected significantly responsive (assessed by testing 
% for a difference of the pooled responses to all stimuli against their baseline
% using a paired Wilcoxon signed-rank test, P<0.05) and selective neurons 
% (significant modulation by one of the stimuli, Kruskal–Wallis test, P<0.05) neurons.
% For the 2,343 neurons that passed both tests, the SNR was calculated as . We observed
% that the SNR distribution was long tailed with a small fraction of cells responding
% with high SNRs. To base the clustering on the clearest signals, we first selected
% the 30% of the cells with largest SNR. Using the Euclidean distance on z scored
% response traces (that is, normalized by their standard deviation), as a similarity
% metric across cells and the ‘furthest distance’ as a measure of distance between
% clusters, we established a hierarchical clustering tree. The tree was thresholded
% to yield 50 different clusters. This method yielded a large number of small clusters,
% which after visual inspection appeared to contain noisy responses (hence very 
% dissimilar to other clusters). We therefore excluded clusters containing <10 cells.
% Applying this criterion, we obtained 13 clusters. Non-clustered cells were then assigned
% to 1 of the 13 clusters with which they had the highest correlation (Pearson correlation
% coefficient) provided that this correlation was >0.1. After this procedure, 1,341 neurons
% were assigned to a cluster while 1,002 cells were not assigned. Inspection of their
% responses showed that the latter were weakly or non-responsive cells.

% %responsiveness
responsiveIdx=cellfun(@(x,y) signrank(x,y), mat2cell(baselines,ones(95,1)),...
    mat2cell(traces,ones(95,1)))<0.05;
% %selectivity . Compare SSD vs NSS? 
% % figure; hold on;  plot(traces(1,:)); plot(baselines(1,:))
% 
% % SNR should be computed beforehand
% % MAD threshold
% % okSNR=max(traces,[],2) > mean(baselines,2)+1.253*mad(baselines,[],2);
% % SNRs=max(traces(responsiveIdx,:),[],2)./mean(baselines(responsiveIdx,:),2);
% % SNRs=max(traces,[],2)./mean(baselines,2);
% % figure; histogram(SNRs,20);
% % % [~,~,foo]= histcounts(SNRs)
% % okSNR=SNRs>1.3;

% zscore traces
zTraces=zscore(traces,[],2);
zTraces=zTraces(responsiveIdx,:); 
zTraces=zTraces(okSNR,:);

figure; plot(zTraces')
figure; plot(traces(responsiveIdx,:)')
figure; plot(traces')
figure; plot(zscore(traces,[],2)')

%Euclidean distance
eDistance = pdist(zTraces); %'cityblock'
% hierarchical clustering tree
clustTree= linkage(eDistance,'complete'); %'average'
% clusters
clusIdx = cluster(clustTree,'maxclust',4); %'criterion','distance','cutoff',2.5*10^4

[proportions,clusterIDs]=hist(clusIdx,unique(clusIdx))

figure;hold on
for clusNum = 1:numel(unique(clusIdx))
    subplot(ceil(numel(unique(clusIdx))/3),3,clusNum); hold on
    clust = clusIdx==clusterIDs(clusNum);
    if size(zTraces(clust,:),1)>1
        plot(nanmean(zTraces(clust,:)))   
    else
        plot((zTraces(clust,:)))   
    end
    yLims=get(gca,'YLim');
%     set(gca,'YLim',[0 max(yLims)]);
        patch([299 300 300 299],[min(yLims) max(yLims) max(yLims) min(yLims)],...
        cmap(clusNum,:),'FaceAlpha',0.5);
    text(1,round(max(yLims)-1),['n=' num2str(size(zTraces(clust,:),1))]);
end

%order by cluster to plot heatmap
[~,clusterSort]=sort(clusIdx);
eDistance = pdist(zTraces(clusterSort,:)); %'cityblock'
figure; imagesc(squareform(eDistance)), colormap(hot); colorbar
figure; heatmap(squareform(eDistance));


%  Distance matrix for the 1,341 clustered neurons.
%  The metric used is d=1?cc, where cc stands for the
%  Pearson correlation coefficient between response traces. 



%% Other method

% Piscopo et al 2013
% http://www.jneurosci.org/content/jneuro/33/11/4642/F2.large.jpg
% This analysis resulted in a set of response parameters for each unit (listed in Table 1).
% Each parameter was zero-centered and normalized to have unity variance across the population,
% and the resulting response profiles were used in a fuzzy k-means clustering algorithm (Matlab Central).
% We heuristically found that setting the fixed number of clusters to n = 6 achieved a tradeoff
% between poor correlations within clusters (too few clusters) and excessive splitting (too many clusters).
% Setting n = 5 merged the “slow” group into the DS/orientation-selective (OS) group and sON/sOFF groups,
% reducing the within-group correlations, and n = 7 split the “slow” group into new groups whose
% characteristic differences were not clear.

% STA RF amplitude	
% Flashed spot polarity	
% Biphasic index	
% Peak flash response latency (ms)	
% Preferred spot diameter (°)	
% Sustain index	
% Correlation of ON/OFF flash resp	
% Preferred speed (°/s)	
% Temporal frequency index	
% Preferred spatial frequency (cpd)	
% Direction selectivity index (DSI)	
% Orientation selectivity index (OSI)	
% Linearity (F1/F0)
% Evoked movie response (sp/s)	
% Spontaneous firing rate (sp/s)	