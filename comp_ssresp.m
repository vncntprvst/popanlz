function [nssResps,nssRespsTrials,sscsResps,sscsRespsTrials,badapl]=comp_ssresp(data,options)
% computes average stop signal response (un-normalized and normalized )
% also outputs individual trial activity for stop signal and baseline epochs

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
%% get stop signal response and baseline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1/ convolve rasters with short_wds ms before stop signal, short_wds after stop signal, 20ms kernel
%time window. Add kernel * 6 ms (see fullgauss_filtconv), e.g. 60 ms at both
% ends, which will be cut.
% data.allndata has 3 column for 3 aligntype. Each cell has 3 or 4 for different conditions
badapl.nsscs=cellfun(@(x) size(x(1).rast,1)==0 | size(x(1).rast,1)==1,data.allndata(:,3));
[nssResps.baseline,nssRespsTrials.baseline]=cellfun(@(x) conv_raster(x(1,1).rast,sigma,...
    0,x(1,1).alignt-(baselineLength+sigma*3),x(1,1).alignt+(sigma*3-1)),...
    data.allndata(~badapl.nsscs,2), 'UniformOutput',false); %500ms period
% no-stop signal trials
[nssResps.short,nssRespsTrials.short]=cellfun(@(x) conv_raster(x(1,1).rast,sigma,...
    0,x(1,1).alignt-(short_wds+sigma*3),x(1,1).alignt+(short_wde+sigma*3)),...
    data.allndata(~badapl.nsscs,3), 'UniformOutput',false); %500ms period
[nssResps.long,nssRespsTrials.long]=cellfun(@(x) conv_raster(x(1,1).rast,sigma,...
    0,x(1,1).alignt-(long_wds+sigma*3),x(1,1).alignt+(short_wde+sigma*3)),...
    data.allndata(~badapl.nsscs,3), 'UniformOutput',false); %700ms period
% stop signal trials with canceled saccades
badapl.sscs=cellfun(@(x) size(x(2).rast,1)==0 | size(x(2).rast,1)==1,data.allndata(:,2)) |...
cellfun(@(x) size(x(2).rast,1)==0 | size(x(2).rast,1)==1,data.allndata(:,3));
[sscsResps.baseline,sscsRespsTrials.baseline]=cellfun(@(x) conv_raster(x(2).rast,sigma,...
    0,x(2).alignt-(baselineLength+sigma*3),x(2).alignt+(sigma*3-1)),...
    data.allndata(~badapl.sscs,2), 'UniformOutput',false); %500ms period
[sscsResps.short,sscsRespsTrials.short]=cellfun(@(x) conv_raster(x(2).rast,sigma,...
    0,x(2).alignt-(short_wds+sigma*3),x(2).alignt+(short_wde+sigma*3)),...
    data.allndata(~badapl.sscs,3), 'UniformOutput',false); %500ms period

% cellfun(@(x) size(x(2).rast,1)>0,data.allndata(:,3));

%% remove bad apples
% badapl=cellfun(@(x) size(x,2)==1, sscsResps.short);
% nssResps.short=nssResps.short(~badapl,:);
nssResps.short=cat(1,nssResps.short{:});
% nssResps.long=nssResps.long(~badapl,:);
nssResps.long=cat(1,nssResps.long{:});
% nssResps.baseline=nssResps.baseline(~badapl,:);
nssResps.baseline=cat(1,nssResps.baseline{:});
% % fullresps=fullresps(~badapl,:);
% % fullresps=cat(1,fullresps{:});
% nssRespsTrials.short=nssRespsTrials.short(~badapl,:);
% nssRespsTrials.long=nssRespsTrials.long(~badapl,:);
% nssRespsTrials.baseline=nssRespsTrials.baseline(~badapl,:);
% 
% sscsResps.short=sscsResps.short(~badapl,:);
sscsResps.short=cat(1,sscsResps.short{:});
% sscsResps.baseline=sscsResps.baseline(~badapl,:);
sscsResps.baseline=cat(1,sscsResps.baseline{:});
% sscsRespsTrials.short=sscsRespsTrials.short(~badapl,:);
% sscsRespsTrials.baseline=sscsRespsTrials.baseline(~badapl,:);

% clusterIdx=clusterIdx(~badapl,:);
% figure; plot(mean(ssresps(clusterIdx==5,:)))

fn = fieldnames(data);
for lp=1:length(fn)
    data.(fn{lp})=data.(fn{lp})(~badapl.nsscs,:);
end

%% normalization
% z-score normalization by baseline - based on pre-target activity
blRespMean=nanmean(nssResps.baseline,2);
blRespSD=nanstd(nssResps.baseline,[],2);
% ssResps.short_blNorm is used for clustering purposes only
nssResps.short_blNorm=(nssResps.short-repmat(blRespMean,1,size(nssResps.short,2)))./repmat(blRespSD,1,size(nssResps.short,2));
nssResps.long_blNorm=(nssResps.long-repmat(blRespMean,1,size(nssResps.long,2)))./repmat(blRespSD,1,size(nssResps.long,2));

% z-score normalization over response period (alternative method, if typically low
% baseline firing rate). Also forces clustering to operate on shapes rather than
% amplitude, by squashing response range
ssRespMean=nanmean(nssResps.short,2);
ssRespSD=nanstd(nssResps.short,[],2);
nssResps.short_respNorm=(nssResps.short-repmat(ssRespMean,1,size(nssResps.short,2)))./repmat(ssRespSD,1,size(nssResps.short,2));
ssRespMean=nanmean(nssResps.long,2);
nssResps.long_respNorm=(nssResps.long-repmat(ssRespMean,1,size(nssResps.long,2)))./repmat(ssRespSD,1,size(nssResps.long,2));

% multitask normalization (see Costello, Salinas, Stanford)
if isfield(data,'normFactor')
    nssResps.short_mtNorm=nssResps.short./repmat(data.normFactor(:,1),1,size(nssResps.short,2));
    nssResps.long_mtNorm=nssResps.long./repmat(data.normFactor(:,1),1,size(nssResps.long,2));
end

%% plot normalized population
% figure;
% subplot(1,3,1)
% imagesc(ssResps.short_blNorm);
% subplot(1,3,2)
% imagesc(ssResps.short_respNorm);
% subplot(1,3,3)
% imagesc(ssResps.short_mtNorm);
% % colormap gray
% % colorbar;
% 
% figure;
% subplot(1,3,1)
% imagesc(ssResps.long_blNorm);
% subplot(1,3,2)
% imagesc(ssResps.long_respNorm);
% subplot(1,3,3)
% imagesc(ssResps.long_mtNorm);
% % colormap parula %copper
% % colorbar;


%% look at "best" cells
% midrange=size(ssResps.short_blNorm,2)/2;
% % Make seed that represent midrange drop / midrange burst / ramp to end / ramp down
% midrangedropseeds=cellfun(@(x) mean(x(1,midrange-150:midrange-50))-mean(x(1,midrange+50:midrange+150)), mat2cell(ssResps.short_blNorm,ones(size(ssResps.short_blNorm,1),1)));
% % ramp to end
% outerrangerampseeds=cellfun(@(x) mean(x(1,length(x)-150:length(x)-1))-mean(x(1,1:150)), mat2cell(ssResps.short_blNorm,ones(size(ssResps.short_blNorm,1),1)));
% % for ramps all the way down, keep only non-bursting / falling response (~monotonic)
% leastdiff_bnorm_ssresps=ssResps.short_blNorm(max(abs(diff(ssResps.short_blNorm)),[],2)<5,:);
% outerrangerampdownseeds=cellfun(@(x) mean(x(1,length(x)-150:length(x)-1))-mean(x(1,1:150)), mat2cell(leastdiff_bnorm_ssresps,ones(size(leastdiff_bnorm_ssresps,1),1)));
% % diff sort works for peaks as well, by opposition, and could be used
% % to separate sharp bursts from smoth bursts (and template 2 from 3 apparently):
% % [~,pkseeds_vals_idx]=sort(max(abs(diff(ssResps.short_blNorm)),[],2),'descend');
% midrangepeakseeds=cellfun(@(x) (mean(x(1,midrange+50:midrange+100))-mean(x(1,midrange-150:midrange-50)))+...
%     (mean(x(1,midrange+50:midrange+100))-mean(x(1,midrange+100:midrange+short_wds))), mat2cell(ssResps.short_blNorm,ones(size(ssResps.short_blNorm,1),1)));
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
