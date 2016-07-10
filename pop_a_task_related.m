function goodrecs=pop_a_task_related(data)
global directory slash;
if isempty(directory)
    [directory,slash]=SetUserDir;
end

%number cells
goodrecs=~cellfun('isempty',data.allsacdelay);

%remove bad apples
fn = fieldnames(data);
for lp=1:length(fn)
    data.(fn{lp})=data.(fn{lp})(goodrecs,:);
end

cellnum=sum(~cellfun('isempty',data.(fn{lp})));

disp([num2str(cellnum) ' cells to analyze'])

%% Convolve traces
sigma=15;
baslineLength=200;
[sacActivity]=cellfun(@(x) conv_raster(x(1,1).rast,sigma,x(1,1).alignt-(200+sigma*3),x(1,1).alignt+(99+sigma*3)), data.allndata(:,1), 'UniformOutput',false); %300ms period
[tgtActivity]=cellfun(@(x) conv_raster(x(1,1).rast,sigma,x(1,1).alignt+(50-sigma*3),x(1,1).alignt+(199+sigma*3)), data.allndata(:,2), 'UniformOutput',false); %250ms period
[bslActivity]=cellfun(@(x) conv_raster(x(1,1).rast,sigma,x(1,1).alignt-(baslineLength+sigma*3),x(1,1).alignt+(sigma*3-1)), data.allndata(:,2), 'UniformOutput',false); %500ms period
% fullresps=cellfun(@(x) conv_raster(x(1,1).rast,sigma,1,size(x(1,1).rast,2)), data.allndata(:,2), 'UniformOutput',false); %full trace

% remove recordings with one trial
properLengthRecs=cellfun(@(x) length(x)>1,bslActivity);
goodrecs(goodrecs)=properLengthRecs;
sacActivity=sacActivity(properLengthRecs);
tgtActivity=tgtActivity(properLengthRecs);
bslActivity=bslActivity(properLengthRecs);

%% Permutation test
% Perform a permutation test on the spike rate in 50 ms intervals
% throughout the defined time period to compare against the baseline
% period.
% For saccade task: 
% visual epoch: 50–200 ms after target onset
% motor epoch:  200 ms before to 100ms after saccade onset
% baseline:     200–150ms prior to target onset.
% Task-related if p value is ?0.01 for two or more of the intervals.

%concatenate test epochs 
testEpochs=cellfun(@(x,y) [x,y], tgtActivity,sacActivity,'UniformOutput', false);

% get mean 50ms baseline 
meanBaseline=cellfun(@(x) mean(reshape(x,[50,4]),2)',bslActivity,'UniformOutput', false);

% run test on each 50ms bins vs mean 50ms baseline
testNum=max(size(testEpochs{1}))/max(size(meanBaseline{1})); 
pValues=cell(1,testNum);
for npttest=1:testNum
    [~, ~, pValues{npttest}] = cellfun(@(x,y) statcond({x(50*(npttest-1)+1:50*(npttest)) y},...
    'method', 'perm', 'naccu', 5000),testEpochs,meanBaseline,'UniformOutput', false);
end

%concatenate p values and get logical index
pValues=[pValues{:}];
sigIdx=sum(cell2mat(pValues)<=0.01,2)>=2;

% map back to initial recording index
goodrecs(goodrecs)=sigIdx;

% [h,p] = ttest(sacActivity{1}(1:50),sacActivity{1}(1:50)+rand(1,50))

% figure;hold on
% plot(sacActivity{1}(1:50))
% plot(sacActivity{1}(1:50)+rand(1,50)-0.5)
% plot(bslActivity{1})
