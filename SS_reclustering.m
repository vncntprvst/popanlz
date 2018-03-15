function [clusidx,clusterIDs]=SS_reclustering(data)
%% old stop-signal re-clustering
% classify according to when peak occurs with respect to stop signal data span
ssd_startstop=[800 700];
clusidx=data.clusters.SaccadeAlignedCluster;
for gsd=1:size(data.allndata,1)
    %         if clusidx(gsd)~=-1
    try
        rasters=data.allndata{gsd,3}(3).rast;
        alignmtt=data.allndata{gsd,3}(3).alignt;
        start=alignmtt-ssd_startstop(1)-half_sixsig; stop=alignmtt+ssd_startstop(2)+half_sixsig;
        normsdf=conv_raster(rasters,conv_sigma,start,stop, data.normFactor(gsd,1)); %normalize by baseline
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
clusterIDs=100+(1:length(unique(clusidx))-1);


% try t-SNE
