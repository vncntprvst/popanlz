function [sacResps,sacRespsTrials,bslRespsTrials,badapl]=comp_sacresp(data,options)
% computes average saccade response (un-normalized and normalized )
% also outputs individual trial activity for saccade and baseline epochs

if exist('options','var')
    sigma=options.sigma;
    baselineLength=options.baselineLength;
    short_wds=options.short_wds;
    short_wde=options.short_wde;
    long_wds=options.long_wds;
else
    sigma=10;
    baselineLength=500;
    short_wds=200;
    short_wde=199;
    long_wds=600;
end
%% get saccade response and baseline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1/ convolve rasters with short_wds ms before saccade, short_wds after saccade, 20ms kernel
%time window. Add kernel * 6 ms (see fullgauss_filtconv), e.g. 60 ms at both
% ends, which will be cut.
% data.allndata has 3 column for 3 aligntype. Each cell has 3 or 4 for different conditions
[bslResps,bslRespsTrials]=cellfun(@(x) conv_raster(x(1,1).rast,sigma,x(1,1).alignt-(baselineLength+sigma*3),x(1,1).alignt+(sigma*3-1)), data.allndata(:,2), 'UniformOutput',false); %500ms period
% fullresps=cellfun(@(x) conv_raster(x(1,1).rast,sigma,1,size(x(1,1).rast,2)), data.allndata(:,2), 'UniformOutput',false); %full response
[sacResps.long,sacRespsTrials.long]=cellfun(@(x) conv_raster(x(1,1).rast,sigma,x(1,1).alignt-(long_wds+sigma*3),x(1,1).alignt+(short_wde+sigma*3)), data.allndata(:,1), 'UniformOutput',false); %400ms period
[sacResps.short,sacRespsTrials.short]=cellfun(@(x) conv_raster(x(1,1).rast,sigma,x(1,1).alignt-(short_wds+sigma*3),x(1,1).alignt+(short_wde+sigma*3)), data.allndata(:,1), 'UniformOutput',false); %400ms period

%% remove bad apples
badapl=cellfun(@(x) size(x,2)==1, sacResps.short);
sacResps.short=sacResps.short(~badapl,:);
sacResps.short=cat(1,sacResps.short{:});
sacResps.long=sacResps.long(~badapl,:);
sacResps.long=cat(1,sacResps.long{:});
bslResps=bslResps(~badapl,:);
bslResps=cat(1,bslResps{:});
% fullresps=fullresps(~badapl,:);
% fullresps=cat(1,fullresps{:});
sacRespsTrials.short=sacRespsTrials.short(~badapl,:);
sacRespsTrials.long=sacRespsTrials.long(~badapl,:);
bslRespsTrials=bslRespsTrials(~badapl,:);
% clusterIdx=clusterIdx(~badapl,:);
% figure; plot(mean(sacresps(clusterIdx==5,:)))

fn = fieldnames(data);
for lp=1:length(fn)
    data.(fn{lp})=data.(fn{lp})(~badapl,:);
end

%% normalization
% z-score normalization by baseline - based on pre-target activity
blRespMean=nanmean(bslResps,2);
blRespSD=nanstd(bslResps,[],2);
% sacResps.short_blNorm is used for clustering purposes only
sacResps.short_blNorm=(sacResps.short-repmat(blRespMean,1,size(sacResps.short,2)))./repmat(blRespSD,1,size(sacResps.short,2));
sacResps.long_blNorm=(sacResps.long-repmat(blRespMean,1,size(sacResps.long,2)))./repmat(blRespSD,1,size(sacResps.long,2));

% z-score normalization over response period (alternative method, if typically low
% baseline firing rate). Also forces clustering to operate on shapes rather than
% amplitude, by squashing response range
sacRespMean=nanmean(sacResps.short,2);
sacRespSD=nanstd(sacResps.short,[],2);
sacResps.short_respNorm=(sacResps.short-repmat(sacRespMean,1,size(sacResps.short,2)))./repmat(sacRespSD,1,size(sacResps.short,2));
sacRespMean=nanmean(sacResps.long,2);
sacResps.long_respNorm=(sacResps.long-repmat(sacRespMean,1,size(sacResps.long,2)))./repmat(sacRespSD,1,size(sacResps.long,2));

% multitask normalization (see Costello, Salinas, Stanford)
if isfield(data,'normFactor')
    sacResps.short_mtNorm=sacResps.short./repmat(data.normFactor(:,1),1,size(sacResps.short,2));
    sacResps.long_mtNorm=sacResps.long./repmat(data.normFactor(:,1),1,size(sacResps.long,2));
end

%% plot normalized population
% figure;
% subplot(1,3,1)
% imagesc(sacResps.short_blNorm);
% subplot(1,3,2)
% imagesc(sacResps.short_respNorm);
% subplot(1,3,3)
% imagesc(sacResps.short_mtNorm);
% % colormap gray
% % colorbar;
% 
% figure;
% subplot(1,3,1)
% imagesc(sacResps.long_blNorm);
% subplot(1,3,2)
% imagesc(sacResps.long_respNorm);
% subplot(1,3,3)
% imagesc(sacResps.long_mtNorm);
% % colormap parula %copper
% % colorbar;


%% look at "best" cells
% midrange=size(sacResps.short_blNorm,2)/2;
% % Make seed that represent midrange drop / midrange burst / ramp to end / ramp down
% midrangedropseeds=cellfun(@(x) mean(x(1,midrange-150:midrange-50))-mean(x(1,midrange+50:midrange+150)), mat2cell(sacResps.short_blNorm,ones(size(sacResps.short_blNorm,1),1)));
% % ramp to end
% outerrangerampseeds=cellfun(@(x) mean(x(1,length(x)-150:length(x)-1))-mean(x(1,1:150)), mat2cell(sacResps.short_blNorm,ones(size(sacResps.short_blNorm,1),1)));
% % for ramps all the way down, keep only non-bursting / falling response (~monotonic)
% leastdiff_bnorm_sacresps=sacResps.short_blNorm(max(abs(diff(sacResps.short_blNorm)),[],2)<5,:);
% outerrangerampdownseeds=cellfun(@(x) mean(x(1,length(x)-150:length(x)-1))-mean(x(1,1:150)), mat2cell(leastdiff_bnorm_sacresps,ones(size(leastdiff_bnorm_sacresps,1),1)));
% % diff sort works for peaks as well, by opposition, and could be used
% % to separate sharp bursts from smoth bursts (and template 2 from 3 apparently):
% % [~,pkseeds_vals_idx]=sort(max(abs(diff(sacResps.short_blNorm)),[],2),'descend');
% midrangepeakseeds=cellfun(@(x) (mean(x(1,midrange+50:midrange+100))-mean(x(1,midrange-150:midrange-50)))+...
%     (mean(x(1,midrange+50:midrange+100))-mean(x(1,midrange+100:midrange+short_wds))), mat2cell(sacResps.short_blNorm,ones(size(sacResps.short_blNorm,1),1)));
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


end
