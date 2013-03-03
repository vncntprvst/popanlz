% Find SSRT for single file
function [overallMeanSSRT,meanIntSSRT,meanSSRT,inhibfun,ssds]=findssrt(recname)
%subjects={'Rigel','Sixx','Hilda'};

% Beta version: best use a recording made with gapstop training and lots of trials
%subject=subjects{subjectnb};
%recname='H53L5A5_20901'; % S148cnttrain
load(recname,'allcodes','alltimes','allbad','saccadeInfo');

%% find stop trials
trialtypes=floor(allcodes(:,2)./10);%./10 only /10 if pooling data together
stoptrials=find(trialtypes==407);
stoptrialcodes=allcodes(stoptrials,:);

%% find canceled and noncanceled stop trials
if find(stoptrialcodes(:,8)==1503,1) %recordings with benchmark code (1503)
    bmssd=1;
    abortedstoptrials=~(floor(stoptrialcodes(:,9)./10)==507);
    noncancel=(stoptrialcodes(~abortedstoptrials,10)==16386) | ... % trials where a saccade is initiated
        (stoptrialcodes(~abortedstoptrials,10)==17385);                % or subsequently breaking fixation
else
    bmssd=0;
    abortedstoptrials=~(floor(stoptrialcodes(:,8)./10)==507);
    noncancel=(stoptrialcodes(~abortedstoptrials,9)==17385) | ...
        (stoptrialcodes(~abortedstoptrials,9)==16386);
end

%% find "desired" delay times (keep in mind that video synch adds ~65ms, so real ssd are variable)

% calculate ssd
stoptrialtimes=alltimes(stoptrials(~abortedstoptrials,:),:);
noncanceltimes=stoptrialtimes(noncancel,:);
canceltimes=stoptrialtimes(~noncancel,:);

if bmssd %benchmark
    nccssd=(noncanceltimes(:,9)-noncanceltimes(:,7))-3; %3 ms added by the state transitions
    ccssd=(canceltimes(:,9)-canceltimes(:,7))-3;
    allssd=nan(size(alltimes,1),1);
    allssd(stoptrials)=(alltimes(stoptrials,9)-alltimes(stoptrials,7))-3;
else
    nccssd=(noncanceltimes(:,8)-noncanceltimes(:,7))-2; %2 ms added by the state transitions
    ccssd=(canceltimes(:,8)-canceltimes(:,7))-2;
end

ssdvalues=unique([ccssd;nccssd]);

%% get saccade delay for signal-respond (noncancel) trial

if bmssd %benchmark
    SRsacdelay=(noncanceltimes(:,10)-noncanceltimes(:,7))-6; %6 ms added by the state transitions
else
    SRsacdelay=(noncanceltimes(:,9)-noncanceltimes(:,7))-6;
end

%% saccade delay for non-stop trials: all good saccade from non-stop trials
%(may yield slightly different results than with left/right parsing method
% used previously)

alllats=reshape({saccadeInfo.latency},size(saccadeInfo));
alllats=alllats';%needs to be transposed because the logical indexing below will be done column by column, not row by row
allgoodsacs=~cellfun('isempty',reshape({saccadeInfo.latency},size(saccadeInfo)));
%removing bad trials
allgoodsacs(logical(allbad),:)=0;
%removing stop trials that may be included
allgoodsacs(floor(allcodes(:,2)./1000)~=6,:)=0;
%indexing good sac trials
% if saccade detection corrected, there may two 'good' saccades
if max(sum(allgoodsacs,2))>1
    twogoods=find(sum(allgoodsacs,2)>1);
    for dblsac=1:length(twogoods)
        allgoodsacs(twogoods(dblsac),find(allgoodsacs(twogoods(dblsac),:),1))=0;
    end
end
sacdelay.all=(cell2mat(alllats(allgoodsacs')))';


%% calculating probabilities for this session
%histo binning method
[~,delaybincenters]=hist([ccssd;nccssd],4);
alldlbincnt=round(delaybincenters');
% sort unique delay bin centers
ssdbins=unique(sort(alldlbincnt));
% narrow ssds to those found more than once
% narssdbins=ssdbins(hist(nccssd,ssdbins)>1); %not for
% individual SSRT calculation
narssdbins=ssdbins(hist(nccssd,ssdbins)>0);
% bin canceled and non-canceled trials according to these ssds
try
    nccssdhist=hist(nccssd,narssdbins);
    ccssdhist=hist(ccssd,narssdbins);
    % find probability to respond
    probaresp=nccssdhist'./(nccssdhist'+ccssdhist');
catch
    probaresp=1;
    narssdbins=0;
end

%diff method
ssdvalues(find(diff(ssdvalues)==1)+1)=ssdvalues(diff(ssdvalues)==1);
ssdvalues=ssdvalues(diff(ssdvalues)>0);
if sum(diff(ssdvalues)==1) % second turn
    ssdvalues(find(diff(ssdvalues)==1))=ssdvalues(diff(ssdvalues)==1)+1;
    ssdvalues=ssdvalues(diff(ssdvalues)>0);
end
% find and keep most prevalent ssds
[ssdtots,ssdtotsidx]=sort((arrayfun(@(x) sum(ccssd==x | ccssd==x-1 | ccssd==x+1),ssdvalues))+...
    (arrayfun(@(x) sum(nccssd==x | nccssd==x-1 | nccssd==x+1),ssdvalues)));
prevssds=sort(ssdvalues(ssdtotsidx(ssdtots>ceil(median(ssdtots))+1)));
try
    if length(prevssds)>=4
        nccssdhist=hist(nccssd,prevssds);
        ccssdhist=hist(ccssd,prevssds);
    else
        nccssdhist=hist(nccssd,ssdvalues);
        ccssdhist=hist(ccssd,ssdvalues);
    end
    % find probability to respond
    probaresp_diff=nccssdhist'./(nccssdhist'+ccssdhist');
    % 	if sum(diff(probaresp_diff)<0) && length(probaresp_diff)-sum(diff(probaresp_diff)<0)>=4
    %         if find(diff(probaresp_diff)>0==0,1)==1
    %             prevssds=prevssds(2:end);
    %             probaresp_diff=probaresp_diff(2:end);
    %         end
    %         prevssds=prevssds(logical([1;diff(probaresp_diff)>0]));
    %         probaresp_diff=probaresp_diff(logical([1;diff(probaresp_diff)>0]));
    %     end
catch
    probaresp_diff=1;
end

%% calculate SSRT
if ~(isempty(narssdbins) || length(narssdbins)==1)
    
    % test monotonicity and keep relevant inhibition function
    if (all(diff(probaresp)>=0) && length(probaresp)>=3) && ~all(diff(probaresp_diff)>=0)
        inhibfun=probaresp;
        ssds=narssdbins;
    elseif all(diff(probaresp_diff)>=0)
        inhibfun=probaresp_diff;
        if length(prevssds)>=4
            ssds=prevssds;
        else
            ssds=ssdvalues;
        end
    else %monotonicity failure
        disp('failed to generate monotonic inhibition function')
        [meanIntSSRT, meanSSRT, overallMeanSSRT, inhibfun, ssds]=deal(NaN);
    end
    
    % calculate SSRT with Boucher et al's method
    try
        [meanIntSSRT, meanSSRT, overallMeanSSRT]= ...
            ssrt_bestfit(sacdelay.all', inhibfun', ssds');
    catch
        sacdelay;
    end
end
end
