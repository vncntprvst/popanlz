function [nssTraces,sscsTraces,ssncsTraces,badapl,blMean]=comp_multiresp(data,options)
% computes average response to different alignments and trial types
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

TestIfNans = @(rasterData,alignment) size(rasterData,1)-...
    sum(isnan(mean(rasterData(:,sigma*3+alignment(1):alignment(2)-3*sigma),2)))<=1;

if ~options.blNorm
    %% no-stop signal trials
    % [nssTraces.baseline]=cellfun(@(x) conv_raster(x(1,1).rast,sigma,...
    %     0,x(1,1).alignt-(long_wds+sigma*3),x(1,1).alignt+(sigma*3-1)),...
    %     {data.target}, 'UniformOutput',false);
    nssTraces.saccade=cellfun(@(x) conv_raster(x(1,1).rast,sigma,...
        0,x(1,1).alignt-(long_wds+sigma*3),x(1,1).alignt+(short_wde+sigma*3)),...
        {data.saccade}, 'UniformOutput',false); % vs sac
    nssTraces.target=cellfun(@(x) conv_raster(x(1,1).rast,sigma,...
        0,x(1,1).alignt-(short_wde+sigma*3),x(1,1).alignt+(long_wds+sigma*3)),...
        {data.target}, 'UniformOutput',false);% vs target
    nssTraces.stopsignal=cellfun(@(x) conv_raster(x(1,1).rast,sigma,...
        0,x(1,1).alignt-(long_wds+sigma*3),x(1,1).alignt+(short_wde+sigma*3)),...
        {data.stopsignal}, 'UniformOutput',false);% vs stop signal
    nssTraces.reward=cellfun(@(x) conv_raster(x(1,1).rast,sigma,...
        0,x(1,1).alignt-(long_wds+sigma*3),x(1,1).alignt+(short_wde+sigma*3)),...
        {data.reward}, 'UniformOutput',false);% vs reward
    
    %% stop signal trials with canceled saccades
    % [sscsTraces.baseline]=cellfun(@(x) conv_raster(x(2).rast,sigma,...
    %     0,x(2).alignt-(long_wds+sigma*3),x(2).alignt+(sigma*3-1)),...
    %     {data.target}, 'UniformOutput',false); % baseline
    sscsTraces.saccade=cellfun(@(x) conv_raster(x(2).rast,sigma,...
        0,x(2).alignt-(long_wds+sigma*3),x(2).alignt+(short_wde+sigma*3)),...
        {data.saccade}, 'UniformOutput',false); % vs sac
    sscsTraces.target=cellfun(@(x) conv_raster(x(2).rast,sigma,...
        0,x(2).alignt-(short_wde+sigma*3),x(2).alignt+(long_wds+sigma*3)),...
        {data.target}, 'UniformOutput',false); % vs target
    sscsTraces.stopsignal=cellfun(@(x) conv_raster(x(3).rast,sigma,...
        0,x(3).alignt-(long_wds+sigma*3),x(3).alignt+(short_wde+sigma*3)),...
        {data.stopsignal}, 'UniformOutput',false); % vs stop signal
    sscsTraces.reward=cellfun(@(x) conv_raster(x(2).rast,sigma,...
        0,x(2).alignt-(long_wds+sigma*3),x(2).alignt+(short_wde+sigma*3)),...
        {data.reward}, 'UniformOutput',false); % vs reward
    
    %% stop signal trials with non-canceled saccades
    badapl=(cellfun(@(rasterData) size({rasterData.rast},2),{data.saccade})<3);

    % [ssncsTraces.baseline]=cellfun(@(x) conv_raster(x(3).rast,sigma,...
    %     0,x(3).alignt-(long_wds+sigma*3),x(3).alignt+(sigma*3-1)),...
    %     {data.target}, 'UniformOutput',false); % baseline
    ssncsTraces.saccade=cellfun(@(x) conv_raster(x(3).rast,sigma,...
        0,x(3).alignt-(long_wds+sigma*3),x(3).alignt+(short_wde+sigma*3)),...
        {data(~badapl).saccade}, 'UniformOutput',false); % vs sac
    ssncsTraces.target=cellfun(@(x) conv_raster(x(3).rast,sigma,...
        0,x(3).alignt-(short_wde+sigma*3),x(3).alignt+(long_wds+sigma*3)),...
        {data(~badapl).target}, 'UniformOutput',false); % vs target
    ssncsTraces.stopsignal=cellfun(@(x) conv_raster(x(4).rast,sigma,...
        0,x(4).alignt-(long_wds+sigma*3),x(4).alignt+(short_wde+sigma*3)),...
        {data(~badapl).stopsignal}, 'UniformOutput',false); % vs stop signal
    ssncsTraces.reward=cellfun(@(x) conv_raster(x(3).rast,sigma,...
        0,x(3).alignt-(long_wds+sigma*3),x(3).alignt+(short_wde+sigma*3)),...
        {data(~badapl).reward}, 'UniformOutput',false); % vs reward
    
    blMean=[];
    
    %% Normalization (baseline substracted)
else
        blMean(:,1)=cellfun(@(x) mean(conv_raster(x(1,1).rast,sigma,...
        0,x(1,1).alignt-(baselineLength+sigma*3),x(1,1).alignt+(sigma*3-1))),...
        {data.target}, 'UniformOutput',false);
        blMean(:,2)=cellfun(@(rasterData) mean(conv_raster(rasterData(2).rast,sigma,...
        0,rasterData(2).alignt-(baselineLength+sigma*3),rasterData(2).alignt+(sigma*3-1))),...
        {data.target},'UniformOutput',false); 
        blMean(:,3)=cellfun(@(rasterData) mean(conv_raster(rasterData(3).rast,sigma,...
        0,rasterData(3).alignt-(baselineLength+sigma*3),rasterData(3).alignt+(sigma*3-1))),...
        {data.target}, 'UniformOutput',false);
        blMean=cellfun(@(x,y,z) nanmean([x,y,z]), blMean(:,1),blMean(:,2),blMean(:,3), 'UniformOutput',false);

%     nssTraces.blMean=cellfun(@(x) mean(conv_raster(x(1,1).rast,sigma,...
%         0,x(1,1).alignt-(long_wds+sigma*3),x(1,1).alignt+(sigma*3-1))),...
%         {data.target}, 'UniformOutput',false);
badapl.nss=cellfun(@(x) size(x(1).rast,1)<=1 || TestIfNans(x(1).rast,[x(1).alignt-...
(short_wds+sigma*3) x(1).alignt+(short_wde+sigma*3)]),{data.stopsignal});

    nssTraces.saccade=cellfun(@(rasterData,blAverage) conv_raster(rasterData(1,1).rast,sigma,...
        0,rasterData(1,1).alignt-(long_wds+sigma*3),rasterData(1,1).alignt+(short_wde+sigma*3))-blAverage,...
        {data.saccade}',blMean, 'UniformOutput',false); % vs sac
    nssTraces.target=cellfun(@(rasterData,blAverage) conv_raster(rasterData(1,1).rast,sigma,...
        0,rasterData(1,1).alignt-(short_wde+sigma*3),rasterData(1,1).alignt+(long_wds+sigma*3))-blAverage,...
        {data.target}',blMean, 'UniformOutput',false);% vs target
    nssTraces.stopsignal=cellfun(@(rasterData,blAverage) conv_raster(rasterData(1,1).rast,sigma,...
        0,rasterData(1,1).alignt-(long_wds+sigma*3),rasterData(1,1).alignt+(short_wde+sigma*3))-blAverage,...
        {data(~badapl.nss).stopsignal}',blMean(~badapl.nss), 'UniformOutput',false);% vs stop signal
    nssTraces.reward=cellfun(@(rasterData,blAverage) conv_raster(rasterData(1,1).rast,sigma,...
        0,rasterData(1,1).alignt-(long_wds+sigma*3),rasterData(1,1).alignt+(short_wde+sigma*3))-blAverage,...
        {data.reward}',blMean, 'UniformOutput',false);% vs reward
    
    %% stop signal trials with canceled saccades
%     sscsTraces.blMean=cellfun(@(rasterData) mean(conv_raster(rasterData(2).rast,sigma,...
%         0,rasterData(2).alignt-(long_wds+sigma*3),rasterData(2).alignt+(sigma*3-1))),...
%         {data.target},'UniformOutput',false); % baseline
badapl.sscs=cellfun(@(x) size(x(2).rast,1)<=1 || TestIfNans(x(2).rast,[x(2).alignt-...
(short_wds+sigma*3) x(2).alignt+(long_wds+sigma*3)]),{data.target}) |...
cellfun(@(x) size(x(3).rast,1)<=1 || TestIfNans(x(3).rast,[x(3).alignt-...
(short_wds+sigma*3) x(3).alignt+(short_wde+sigma*3)]),{data.stopsignal});

    sscsTraces.saccade=cellfun(@(rasterData,blAverage) conv_raster(rasterData(2).rast,sigma,...
        0,rasterData(2).alignt-(long_wds+sigma*3),rasterData(2).alignt+(short_wde+sigma*3))-blAverage,...
        {data.saccade}',blMean, 'UniformOutput',false); % vs sac
    sscsTraces.target=cellfun(@(rasterData,blAverage) conv_raster(rasterData(2).rast,sigma,...
        0,rasterData(2).alignt-(short_wde+sigma*3),rasterData(2).alignt+(long_wds+sigma*3))-blAverage,...
        {data.target}',blMean, 'UniformOutput',false); % vs target
    sscsTraces.stopsignal=cellfun(@(rasterData,blAverage) conv_raster(rasterData(3).rast,sigma,...
        0,rasterData(3).alignt-(long_wds+sigma*3),rasterData(3).alignt+(short_wde+sigma*3))-blAverage,...
        {data.stopsignal}',blMean, 'UniformOutput',false); % vs stop signal
    sscsTraces.reward=cellfun(@(rasterData,blAverage) conv_raster(rasterData(2).rast,sigma,...
        0,rasterData(2).alignt-(long_wds+sigma*3),rasterData(2).alignt+(short_wde+sigma*3))-blAverage,...
        {data.reward}',blMean, 'UniformOutput',false); % vs reward
    
    %% stop signal trials with non-canceled saccades
    badapl.ssncs=cellfun(@(x) size(x(3).rast,1)<=1 || TestIfNans(x(3).rast,[x(3).alignt-...
(short_wds+sigma*3) x(3).alignt+(long_wds+sigma*3)]),{data.target}) |...
cellfun(@(x) size(x(4).rast,1)<=1 || TestIfNans(x(4).rast,[x(4).alignt-...
(short_wds+sigma*3) x(4).alignt+(short_wde+sigma*3)]),{data.stopsignal});
%     badapl.ssncs=(cellfun(@(rasterData) size({rasterData.rast},2),{data.saccade})<3);
%     ssncsTraces.blMean=cellfun(@(rasterData) mean(conv_raster(rasterData(3).rast,sigma,...
%         0,rasterData(3).alignt-(long_wds+sigma*3),rasterData(3).alignt+(sigma*3-1))),...
%         {data.target}, 'UniformOutput',false); % baseline
    ssncsTraces.saccade=cellfun(@(rasterData,blAverage) conv_raster(rasterData(3).rast,sigma,...
        0,rasterData(3).alignt-(long_wds+sigma*3),rasterData(3).alignt+(short_wde+sigma*3))-blAverage,...
        {data(~badapl.ssncs).saccade}',blMean(~badapl.ssncs), 'UniformOutput',false); % vs sac
    ssncsTraces.target=cellfun(@(rasterData,blAverage) conv_raster(rasterData(3).rast,sigma,...
        0,rasterData(3).alignt-(short_wde+sigma*3),rasterData(3).alignt+(long_wds+sigma*3))-blAverage,...
        {data.target}',blMean, 'UniformOutput',false); % vs target
    ssncsTraces.stopsignal=cellfun(@(rasterData,blAverage) conv_raster(rasterData(4).rast,sigma,...
        0,rasterData(4).alignt-(long_wds+sigma*3),rasterData(4).alignt+(short_wde+sigma*3))-blAverage,...
        {data.stopsignal}',blMean, 'UniformOutput',false); % vs stop signal
    ssncsTraces.reward=cellfun(@(rasterData,blAverage) conv_raster(rasterData(3).rast,sigma,...
        0,rasterData(3).alignt-(long_wds+sigma*3),rasterData(3).alignt+(short_wde+sigma*3))-blAverage,...
        {data.reward}',blMean, 'UniformOutput',false); % vs reward
    
    
    %     fieldNames={'saccade';'target';'stopsignal';'reward'};
    %     baseline=vertcat(nssTraces.baseline{:}); %,...
    %     %         sscsTraces.baseline{:}, ssncsTraces.baseline{:});
    %     blRespMean=nanmean(baseline,2);
    %     % blRespSD=nanstd(nssResps.baseline,[],2);
    %     for fn=[1,2,4]
    %         (fieldNames{fn})=mat2cell(vertcat(...
    %             nssTraces.(fieldNames{fn}){:})-blRespMean,ones(100,1));
    %     end
    %     for fn=[2,3,4]
    %         sscsTraces.(fieldNames{fn})=mat2cell(...
    %             vertcat(sscsTraces.(fieldNames{fn}){:})-blRespMean,ones(100,1));
    %     end
    %     for fn=[1,2,3]
    %         ssncsTraces.(fieldNames{fn})=mat2cell(...
    %             vertcat(ssncsTraces.(fieldNames{fn}){:})-blRespMean(~badapl),ones(99,1));
    %     end
end

%% cells with not enough data in one of the conditions
badapl.all = badapl.nss | badapl.sscs | badapl.ssncs;
% cellfun(@(x,y,z,u,v,w,s,t) size(x,2)==1 || size(y,2)==1 || size(z,2)==1 ||...
%     size(u,2)==1 || size(v,2)==1 || size(w,2)==1 || size(s,2)==1 || size(t,2)==1, ...
%     nssTraces.saccade, nssTraces.target, nssTraces.reward, sscsTraces.target,...
%     sscsTraces.stopsignal, sscsTraces.reward, nssTraces.stopsignal, sscsTraces.saccade)...
%     | badapl;