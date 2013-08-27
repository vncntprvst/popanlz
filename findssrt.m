% Find SSRT for single file
function [overallMeanSSRT,meanIntSSRT,meanSSRT,inhibfun,ssds,tachomc,tachowidth]=findssrt(recname, plots)
global directory;

if nargin < 2 | isempty(plots)
     plots = 0;
end

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

        %keep track of trial number
        goodstoptrials=stoptrials(~abortedstoptrials);
        noncanceltrials=goodstoptrials(noncancel);
        canceltrials=goodstoptrials(~noncancel);
        allsactrials=find(trialtypes==604);
        if find(stoptrialcodes(:,8)==1503,1)
            sactrials=allsactrials(floor(allcodes(allsactrials,9)./10)==704);
        else
            sactrials=allsactrials(floor(allcodes(allsactrials,8)./10)==704);
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
%keeping sac info of non-canceled SS trials
        allncsacs=allgoodsacs;
        allncsacs(floor(allcodes(:,2)./1000)==6,:)=0; % nullifying NSS trials
        nasstrials=stoptrials(~abortedstoptrials);
        allncsacs(nasstrials(~noncancel),:)=0; % nullifying CSS trials
        if max(sum(allncsacs,2))>1
            twogoods=find(sum(allncsacs,2)>1);
            for dblsac=1:length(twogoods)
                allncsacs(twogoods(dblsac),find(allncsacs(twogoods(dblsac),:),1))=0;
            end
        end
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
ncsacdelay=cell2mat(alllats(allncsacs'));

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

 %% get tachometric curve. Build matrix with 1. SSDs 2. RTs 3. success
        SSDRTs=inf(length(allbad),3); % may differ from sum(~allbad)
        ccssd(ismember(ccssd,unique(ccssd-1)))=ccssd(ismember(ccssd,unique(ccssd-1)))+1;
        SSDRTs(canceltrials,1)=ccssd;
        nccssd(ismember(nccssd,unique(nccssd-1)))=nccssd(ismember(nccssd,unique(nccssd-1)))+1;
        SSDRTs(logical(sum(allncsacs,2)),1)=nccssd(ismember(goodstoptrials(noncancel),find(sum(allncsacs,2)))); %got to remove some ssd when sacs are not available
        SSDRTs(logical(sum(allgoodsacs,2)),2)=sacdelay.all;
        SSDRTs(logical(sum(allncsacs,2)),2)=ncsacdelay;
        SSDRTs(logical(sum(allncsacs,2)),3)=0;
        SSDRTs(logical(allbad),3)=0;
        SSDRTs([find(sum(allgoodsacs,2));canceltrials],3)=1;
        try
            [tachomc xtach tach rPTc rPTe] = tachCM2(SSDRTs);
        catch
            [tachomc xtach tach rPTc rPTe] = deal(NaN);
        end
        
        if plots
            tachoh=figure;
            filttach=gauss_filtconv(tach,4);
            coretach=filttach(xtach>-20 & xtach<170);
            tachowidth=length(coretach(coretach>=0.1 & coretach<=0.9));
            plot(xtach(xtach>-20 & xtach<170),coretach,'LineWidth',3);
            title('Tachometric curve','FontName','calibri','FontSize',15);
            hxlabel=xlabel(gca,'rPT (ms)','FontName','calibri','FontSize',12);
            set(gca,'Xlim',[-20 170],'XTick',[0:50:150],'TickDir','out','box','off'); %'XTickLabel',[50:50:400]
            hylabel=ylabel(gca,'Fraction cancelled','FontName','calibri','FontSize',12);
            set(gca,'Ylim',[0 1],'TickDir','out','box','off');
            exportfigname=[directory,'figures\cmd\',recname,'_tacho'];
            plot2svg([exportfigname,'.svg'],gcf, 'png');
            delete(tachoh);
            
            sPTdistribh=figure;
            filtrPTc = gauss_filtconv(rPTc,10);
            filtrPTe = gauss_filtconv(rPTe,10);
            plot(xtach,filtrPTc,'r','LineWidth',3);            
            hold on
%           plot(xtach,rPTc,'r')
%           plot(xtach,rPTe,'b')
            plot(xtach,filtrPTe,'b','LineWidth',3);
            title('Distribution of rPTs','FontName','calibri','FontSize',15);
            hxlabel=xlabel(gca,'rPT (ms)','FontName','calibri','FontSize',12);
            set(gca,'Xlim',[min(xtach) max(xtach)],'XTick',[min(xtach):50:max(xtach)],'TickDir','out','box','off'); %'XTickLabel',[50:50:400]
            hylabel=ylabel(gca,'Normalized frequency','FontName','calibri','FontSize',12);
            set(gca,'TickDir','out','box','off');
            set(gca,'XTick',[-200:100:500]);
            set(gca,'XTickLabel',[-200:100:500]);
            legend('rPTc','rPTe');
            exportfigname=[directory,'figures\cmd\',recname,'_rPT'];
            plot2svg([exportfigname,'.svg'],gcf, 'png');
            delete(sPTdistribh);
        else
            tachowidth=NaN;
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
else
    [meanIntSSRT, meanSSRT, overallMeanSSRT, inhibfun, ssds]=deal(NaN);
end
end
