function [sacresps,bnorm_sacresps,rnorm_sacresps,sacrespsTrials,bslrespsTrials]=comp_sacresp(data,dataField)
% computes average saccade response (un-normalized and normalized )
% also outputs individual trial activity for saccade and baseline epochs

%% get saccade response and baseline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1/ convolve rasters with 200ms before saccade, 200 after saccade, 20ms kernel
%time window. Add kernel * 6 ms (see fullgauss_filtconv), e.g. 60 ms at both
% ends, which will be cut.
% data.(dataField).allndata has 3 column for 3 aligntype. Each cell has 3 or 4 for different conditions
sigma=10;
baslineLength=500;
[sacresps,sacrespsTrials]=cellfun(@(x) conv_raster(x(1,1).rast,sigma,x(1,1).alignt-(200+sigma*3),x(1,1).alignt+(199+sigma*3)), data.(dataField).allndata(:,1), 'UniformOutput',false); %400ms period
[bslresps,bslrespsTrials]=cellfun(@(x) conv_raster(x(1,1).rast,sigma,x(1,1).alignt-(baslineLength+sigma*3),x(1,1).alignt+(sigma*3-1)), data.(dataField).allndata(:,2), 'UniformOutput',false); %500ms period
fullresps=cellfun(@(x) conv_raster(x(1,1).rast,sigma,1,size(x(1,1).rast,2)), data.(dataField).allndata(:,2), 'UniformOutput',false); %full response
%% remove bad apples
badapl=cellfun(@(x) size(x,2)==1, sacresps);
sacresps=sacresps(~badapl,:);
sacresps=cat(1,sacresps{:});
bslresps=bslresps(~badapl,:);
bslresps=cat(1,bslresps{:});
fullresps=fullresps(~badapl,:);
% fullresps=cat(1,fullresps{:});
sacrespsTrials=sacrespsTrials(~badapl,:);
bslrespsTrials=bslrespsTrials(~badapl,:);
% clusterIdx=clusterIdx(~badapl,:);
% figure; plot(mean(sacresps(clusterIdx==5,:)))

data.(dataField).allalignmnt=data.(dataField).allalignmnt(~badapl,:);
data.(dataField).allmssrt_tacho=data.(dataField).allmssrt_tacho(~badapl,1);
%allpk=allpk(~badapl,:); %not needed
data.(dataField).allndata=data.(dataField).allndata(~badapl,:);
%all_rec_id=all_rec_id(~badapl,1); %not needed
%allstats=allstats(~badapl,1); %not needed
data.(dataField).allprevssd=data.(dataField).allprevssd(~badapl,:);
data.(dataField).allssds=data.(dataField).allssds(~badapl,:);
data.(dataField).allsacdelay=data.(dataField).allsacdelay(~badapl,:);
data.(dataField).allprefdir=data.(dataField).allprefdir(~badapl,:);
%alltrialidx=alltrialidx(~badapl,:); %not needed
data.(dataField).alldb=data.(dataField).alldb(~badapl,:);
%% normalization
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
%     align=data.(dataField).allndata{top_drop(topfig), 3}(4).alignt;
%     rasters=((data.(dataField).allndata{top_drop(topfig), 3}(4).rast(:,align-800:align+800)));
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
