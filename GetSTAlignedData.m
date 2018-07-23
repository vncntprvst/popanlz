function alignedData=GetSTAlignedData(fileName, aligntype, globVars,allcodes)
global rexnumtrials recDir
rexnumtrials = globVars{1};
recDir = globVars{2};
% recDir='E:\Data\Recordings\processed\Sixx\';

tasktype='st_saccades';

%get trial directions
trialtypes=allcodes(:,2);
trialdirs=unique(trialtypes-floor(trialtypes./10)*10);

%alignments=1:3;
switch aligntype
    case 'saccade' %sac vs stop
        firstalign=6;
        secondalign=4;
        option=NaN;
    case 'target'% tgt vs stop
        firstalign=7;
        secondalign=[];
        option=NaN;
    case 'reward' % reward
        firstalign=4;
        secondalign=[];
        option=NaN;
end

includebad=0;
spikechannel=1; %select appropriate cluster
keepdir='compall'; %alldir %which sac directions
togrey=[];
singlerastplot=0;

%% use GUI-independent prealign
alignedData = prealign(fileName, trialdirs, tasktype, firstalign,...
    secondalign,  includebad, spikechannel, keepdir,...
    togrey, singlerastplot, option);