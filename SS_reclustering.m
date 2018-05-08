 function [clusidx,clusterIDs]=SS_reclustering(dataset)
%% old stop-signal re-clustering
% classify according to when peak occurs with respect to stop signal data span
ssd_startstop=[800 700];
% dataset=data.gsdata;
clusidx=dataset.clusters.SaccadeAlignedCluster;
for gsd=1:size(dataset.allndata,1)
    %         if clusidx(gsd)~=-1
    try
        rasters=dataset.allndata{gsd,3}(3).rast;
        alignmtt=dataset.allndata{gsd,3}(3).alignt;
        start=alignmtt-ssd_startstop(1)-half_sixsig; stop=alignmtt+ssd_startstop(2)+half_sixsig;
        normsdf=conv_raster(rasters,conv_sigma,start,stop, dataset.normFactor(gsd,1)); %normalize by baseline
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
    %         end
end
% clusterIDs=100+(1:length(unique(clusidx))-1);

clusidx(clusidx==-1)=102;
clusidx(clusidx==2)=1;
clusidx(clusidx==101)=2;
clusidx(clusidx==102)=3;
clusidx(clusidx==103)=4;

if ~isfield(dataset,'long') 
% get traces
conv_sigma=50;
half_sixsig=conv_sigma*3; %half kernel window
ssd_startstop=[500 300];
clusidx=dataset.clusters.SaccadeAlignedCluster;
normsdf=nan(size(dataset.allndata,1),sum(ssd_startstop)+1);
for gsd=1:size(dataset.allndata,1)
    try
        rasters=dataset.allndata{gsd,3}(3).rast;
        alignmtt=dataset.allndata{gsd,3}(3).alignt;
        start=alignmtt-ssd_startstop(1)-half_sixsig; stop=alignmtt+ssd_startstop(2)+half_sixsig;
        normsdf(gsd,:)=conv_raster(rasters,conv_sigma,0,start,stop); %,dataset.normFactor(gsd,1));
    catch
        continue
    end
end
badApples=isnan(sum(normsdf,2));
labels=dataset.clusters.Profile;
labels=labels(~badApples);
else 
    normsdf=dataset.long; %long_blNorm;
end

cmap=lines; cmap=[0,0,0;cmap];

% comps = tsne(normsdf,'NumPCAComponents',5,'Algorithm','exact','Distance','cosine');
% figure; gscatter(comps(:,1), comps(:,2), labels ,cmap,'.',7);

[~,score] = pca(single(normsdf));
% figure; scatter3(score(:,1),score(:,2),score(:,3))

cbD = pdist(normsdf,'cityblock'); %'cityblock'
clustTreeCityBlock = linkage(cbD,'average'); %'average'
% figure;[h,nodes] = dendrogram(clustTreeCityBlock,12);

cosD = pdist(score(:,1:4),'cosine');
clustTreeCos = linkage(cosD,'average');
% figure;[h,nodes] = dendrogram(clustTreeCos,12);

%cluster
clusidx = cluster(clustTreeCityBlock,'criterion','distance','cutoff',2.5*10^4);

clusidx = cluster(clustTreeCos,'criterion','distance','cutoff',.5);
[proportions,clusterIDs]=hist(clusidx,unique(clusidx))

figure;hold on
for clusNum = 1:numel(unique(clusidx))
    subplot(ceil(numel(unique(clusidx))/3),3,clusNum); hold on
    clust = clusidx==clusterIDs(clusNum);
    if size(normsdf(clust,:),1)>1
        plot(nanmean(normsdf(clust,:)))   
    else
        plot((normsdf(clust,:)))   
    end
    yLims=get(gca,'YLim');
%     set(gca,'YLim',[0 max(yLims)]);
        patch([499 500 500 499],[min(yLims) max(yLims) max(yLims) min(yLims)],...
        cmap(clusNum,:),'FaceAlpha',0.5);
    text(1,round(max(yLims)-5),['n=' num2str(size(normsdf(clust,:),1))]);
end
title('HC on Traces, no norm, Cityblock')

data.gsdata.clusters = [data.gsdata.clusters  array2table(clusidx,'VariableNames',{'StopSignalAlignedCluster'})];

