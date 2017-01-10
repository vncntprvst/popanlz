% Find SSRT for single file
function [mssrt,inhibfun,ccssd,nccssd,ssds,...
    tachomc,tachowidth,sacdelay,rewtimes,prevssds,trialidx,alltacho]=findssrt(recname, plots)
global directory;

if nargin < 2 | isempty(plots)
    plots = 0;
end

% Beta version: best use a recording made with gapstop training and lots of trials
%subject=subjects{subjectnb};
%recname='H53L5A5_20901'; % S148cnttrain
load(recname,'allcodes','alltimes','allbad','saccadeInfo');
if ~isfield(saccadeInfo,'latency')
    [mssrt,inhibfun,ccssd,nccssd,ssds,tachomc,tachowidth,rewtimes,prevssds,trialidx,...
        sacdelay.nsst,sacdelay.ncst,sacdelay.encst]=deal(NaN);
    return;
end
%% find stop trials
trialtypes=floor(allcodes(:,2)./10);%./10 only /10 if pooling data together
trialidx.stoptrials=find(trialtypes==407);
stoptrialcodes=allcodes(trialidx.stoptrials,:);

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

%% saccade delay for non-stop trials: all good saccade from non-stop trials
%(may yield slightly different results than with left/right parsing method
% used previously)

alllats=reshape({saccadeInfo.latency},size(saccadeInfo));
alllats=alllats';%needs to be transposed because the logical indexing below will be done column by column, not row by row
allgoodsacs=~cellfun('isempty',reshape({saccadeInfo.latency},size(saccadeInfo)));
%weeding out bad trials that are not stop trials
allgoodsacs(logical(allbad)'&trialtypes~=407,:)=0;
%also removing those weird trial with 2222 code
[~,badcol]=ind2sub(size(allcodes),find(allcodes==2222,1));
allgoodsacs(allcodes(:,badcol)==2222,:)=0;
%keeping sac info of non-canceled SS trials
allncsacs=allgoodsacs;
allncsacs(floor(allcodes(:,2)./1000)==6,:)=0; % nullifying NSS trials
nasstrials=trialidx.stoptrials(~abortedstoptrials);
allncsacs(nasstrials(~noncancel),:)=0; % nullifying CSS trials
if max(sum(allncsacs,2))>1
    twogoods=find(sum(allncsacs,2)>1);
    for dblsac=1:length(twogoods)
        allncsacs(twogoods(dblsac),find(allncsacs(twogoods(dblsac),:),1))=0;
    end
end
%removing stop trials that may be included
allgoodsacs(trialtypes==407,:)=0;
%indexing good sac trials
% if saccade detection corrected, there may two 'good' saccades
if max(sum(allgoodsacs,2))>1
    twogoods=find(sum(allgoodsacs,2)>1);
    for dblsac=1:length(twogoods)
        allgoodsacs(twogoods(dblsac),find(allgoodsacs(twogoods(dblsac),:),1))=0;
    end
end
sacdelay.nsst=(cell2mat(alllats(allgoodsacs')))';
%get reward time for NSS trials
goodsactimes=alltimes(logical(sum(allgoodsacs,2)),:);
trialidx.goodsac=find(sum(allgoodsacs,2));
rewtimes=goodsactimes(allcodes(logical(sum(allgoodsacs,2)),:)==1030);

%% find trials to the few early saccades that happened just before SSD
aborsstrials=trialidx.stoptrials(abortedstoptrials);
earlyncsac=aborsstrials(ismember(aborsstrials,find(sum(allncsacs,2))));

% To adjust abortedstoptrials and noncancel (but impractical because can't
% calculate ssd from those early trials:
% if ~isempty(earlyncsac)
%     abortedstoptrials(ismember(trialidx.stoptrials,earlyncsac))=0;
%     if find(stoptrialcodes(:,8)==1503,1)
%         noncancel=logical(sum((stoptrialcodes(~abortedstoptrials,9:10)==17385) |...
%         (stoptrialcodes(~abortedstoptrials,9:10)==16386),2));
%     else
%         noncancel=logical(sum((stoptrialcodes(~abortedstoptrials,8:9)==17385) |...
%         (stoptrialcodes(~abortedstoptrials,8:9)==16386),2));
%     end
% end

%Instead, adjust allncsacs
sacdelay.encst=([alllats{:,earlyncsac}]);
if length(sacdelay.encst)~=length(earlyncsac)
    dbstack
    return
end
allncsacs(earlyncsac,:)=0;
sacdelay.ncst=cell2mat(alllats(allncsacs'));

%% keep track of trial number
goodstoptrials=trialidx.stoptrials(~abortedstoptrials);
trialidx.noncancel=goodstoptrials(noncancel);
trialidx.cancel=goodstoptrials(~noncancel);
% allsactrials=find(trialtypes==604);
% if find(stoptrialcodes(:,8)==1503,1)
%     sactrials=allsactrials(floor(allcodes(allsactrials,9)./10)==704);
%     if isempty(sactrials) && find(floor(allcodes(allsactrials,8)./10)==704,1) %rare case where 1503 was not added everywhere
%         sactrials=allsactrials(floor(allcodes(allsactrials,8)./10)==704);
%     end
% else
%     sactrials=allsactrials(floor(allcodes(allsactrials,8)./10)==704);
% end

%% find "desired" delay times (keep in mind that video synch adds ~65ms, so real ssd are variable)

% calculate ssd
stoptrialtimes=alltimes(trialidx.stoptrials(~abortedstoptrials,:),:);
noncanceltimes=stoptrialtimes(noncancel,:);
canceltimes=stoptrialtimes(~noncancel,:);

if bmssd %benchmark
    nccssd=(noncanceltimes(:,9)-noncanceltimes(:,7))-3; %3 ms added by the state transitions
    ccssd=(canceltimes(:,9)-canceltimes(:,7))-3;
    allssd=(alltimes(trialidx.stoptrials(~abortedstoptrials),9)-alltimes(trialidx.stoptrials(~abortedstoptrials),7))-3;
else
    nccssd=(noncanceltimes(:,8)-noncanceltimes(:,7))-2; %2 ms added by the state transitions
    ccssd=(canceltimes(:,8)-canceltimes(:,7))-2;
    allssd=(alltimes(trialidx.stoptrials(~abortedstoptrials),8)-alltimes(trialidx.stoptrials(~abortedstoptrials),7))-3;
end

ssdvalues=unique([ccssd;nccssd]);

%% get saccade delay for signal-respond (noncancel) trial

if bmssd %benchmark
    SRsacdelay=(noncanceltimes(:,10)-noncanceltimes(:,7))-6; %6 ms added by the state transitions
else
    SRsacdelay=(noncanceltimes(:,9)-noncanceltimes(:,7))-6;
end



%% calculating probabilities for this session
%histo binning method: thing of the past
% [~,delaybincenters]=hist([ccssd;nccssd],4);
% alldlbincnt=round(delaybincenters');
% % sort unique delay bin centers
% ssdbins=unique(sort(alldlbincnt));
% % narrow ssds to those found more than once
% % narssdbins=ssdbins(hist(nccssd,ssdbins)>1); %not for
% % individual SSRT calculation
% narssdbins=ssdbins(hist(nccssd,ssdbins)>0);
% % bin canceled and non-canceled trials according to these ssds
% try
%     nccssdhist=hist(nccssd,narssdbins);
%     ccssdhist=hist(ccssd,narssdbins);
%     % find probability to respond
%     probaresp=nccssdhist'./(nccssdhist'+ccssdhist');
% catch
%     probaresp=1;
%     narssdbins=0;
% end

%diff method
% ssdvalues(find(diff(ssdvalues)==1)+1)=ssdvalues(diff(ssdvalues)==1);
% ssdvalues=ssdvalues(diff(ssdvalues)>0);
% if sum(diff(ssdvalues)==1) % second turn
%     ssdvalues(find(diff(ssdvalues)==1))=ssdvalues(diff(ssdvalues)==1)+1;
%     ssdvalues=ssdvalues(diff(ssdvalues)>0);
% end

while sum(diff(ssdvalues)==1)
    ssdvalues(diff(ssdvalues)==1)=ssdvalues(diff(ssdvalues)==1)+1;
end
ssdvalues=unique(ssdvalues);

%% find and keep most prevalent ssds
[ssdtots,ssdtotsidx]=sort((arrayfun(@(x) sum(ccssd<=x+3 & ccssd>=x-3),ssdvalues))); %+...
%     (arrayfun(@(x) sum(nccssd<=x+3 & nccssd>=x-3),ssdvalues)));
if length(ssdtots)>3 && max(ssdtots)<4 && sum(diff(ssdtots(end-3:end)))==1
    prevssds=sort(ssdvalues(ssdtotsidx(length(ssdtotsidx)-2:length(ssdtotsidx))));
else
    prevssds=sort(ssdvalues(ssdtotsidx(ssdtots>(median(ssdtots(ssdtots>1))-std(ssdtots)))));
end
raise=1;
while length(prevssds)>4 && raise<4
    prevssds=sort(ssdvalues(ssdtotsidx(ssdtots>(median(ssdtots)-std(ssdtots)+raise))));
    raise=raise+1;
end

%% cancellation probability (aka inhibition function)
try
    if length(prevssds)>=2
        nccssdhist=hist(nccssd,sort(prevssds));
        ccssdhist=hist(ccssd,sort(prevssds));
        if sum(diff(ccssdhist)<0)>=length(ccssdhist)-2 && diff(nccssdhist(end-1:end))<0
            prevssds=prevssds(1:end-1);
            nccssdhist=hist(nccssd,sort(prevssds));
            ccssdhist=hist(ccssd,sort(prevssds));
        end
        nccssdhist=nccssdhist(nccssdhist>0);
        ccssdhist=ccssdhist(nccssdhist>0);
        prevssds=prevssds(nccssdhist>0);
    else
        nccssdhist=hist(nccssd,sort(ssdvalues(ssdtotsidx(ssdtots>0))));
        ccssdhist=hist(ccssd,sort(ssdvalues(ssdtotsidx(ssdtots>0))));
        nccssdhist=nccssdhist(nccssdhist>0);
        ccssdhist=ccssdhist(nccssdhist>0);
        ssdvalues=ssdvalues(nccssdhist>0);
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
SSDRTs(trialidx.cancel,1)=ccssd;
nccssd(ismember(nccssd,unique(nccssd-1)))=nccssd(ismember(nccssd,unique(nccssd-1)))+1;
SSDRTs(logical(sum(allncsacs,2)),1)=nccssd(ismember(goodstoptrials(noncancel),find(sum(allncsacs,2)))); %got to remove some ssd when sacs are not available
SSDRTs(logical(sum(allgoodsacs,2)),2)=sacdelay.nsst;
SSDRTs(logical(sum(allncsacs,2)),2)=sacdelay.ncst;
SSDRTs(logical(sum(allncsacs,2)),3)=0;
SSDRTs(logical(allbad),3)=0;
SSDRTs([find(sum(allgoodsacs,2));trialidx.cancel],3)=1;
try
    [tachomc, xtach, tach, rPTc, rPTe] = tachCM2(SSDRTs);
catch
    [tachomc, xtach, tach, rPTc, rPTe] = deal(NaN);
end
if ~isnan(xtach)
    if sum(tach(find(tach==1,1):end)<1)
%         close all; figure; hold on; plot(xtach,tach)
        tach(end:-1:find(tach==1,1))=1;
%         plot(xtach,tach); hold off
    end
end
    % smooth tacho curve and get tachowidth
filttach=gauss_filtconv(tach,4);

coretach=filttach(xtach>-20 & xtach<170);
tachowidth=length(coretach(coretach>=0.1 & coretach<=0.9));
if xtach(find(tach==max(tach),1))>50 & xtach(find(tach==max(tach),1))<120
    alltacho.tachomc=tachomc;
    alltacho.tachowidth=tachowidth;
    alltacho.xtach=xtach;
    alltacho.tach=tach;
    alltacho.rPTc=rPTc;
    alltacho.rPTe=rPTe;
else
	[alltacho.tachomc,alltacho.tachowidth,alltacho.xtach,...
        alltacho.tach,alltacho.rPTc,alltacho.rPTe]=deal(NaN);
end

if plots
    psychoplots=figure('color','white','position',[2315	200	524	636]);
    subplot(3,2,1)
    plot(xtach(xtach>-20 & xtach<170),coretach,'LineWidth',2);
    title('Tachometric curve','FontName','calibri','FontSize',12);
    hxlabel=xlabel(gca,'rPT (ms)','FontName','calibri','FontSize',12);
    set(gca,'Xlim',[-20 170],'XTick',[0:50:150],'TickDir','out','box','off'); %'XTickLabel',[50:50:400]
    hylabel=ylabel(gca,'Fraction cancelled','FontName','calibri','FontSize',12);
    set(gca,'Ylim',[0 1],'TickDir','out','box','off');
    
    subplot(3,2,2)
    filtrPTc = gauss_filtconv(rPTc,10);
    filtrPTe = gauss_filtconv(rPTe,10);
    plot(xtach,filtrPTc,'r','LineWidth',2);
    hold on
    %           plot(xtach,rPTc,'r')
    %           plot(xtach,rPTe,'b')
    plot(xtach,filtrPTe,'b','LineWidth',2);
    title('Distribution of rPTs','FontName','calibri','FontSize',15);
    hxlabel=xlabel(gca,'rPT (ms)','FontName','calibri','FontSize',12);
    set(gca,'Xlim',[min(xtach) max(xtach)],'XTick',[min(xtach):50:max(xtach)],'TickDir','out','box','off'); %'XTickLabel',[50:50:400]
    hylabel=ylabel(gca,'Normalized frequency','FontName','calibri','FontSize',12);
    set(gca,'TickDir','out','box','off');
    set(gca,'XTick',[-200:100:500]);
    set(gca,'XTickLabel',[-200:100:500]);
    %             legend('rPTc','rPTe');
    
    subplot(3,2,3:6)
    delaydistribedges=min(sacdelay.nsst)-1:10:max(sacdelay.nsst)+1;
    delaydistribhhist=histc(sort(sacdelay.nsst),delaydistribedges);
    delayfreq=delaydistribhhist./sum(delaydistribhhist);
    delayfreq=fullgauss_filtconv(delayfreq,2,0);
    plot(delaydistribedges,delayfreq,'-k','LineWidth',2);
    xlabel('Saccade delay (ms)','FontName','calibri','FontSize',12);
    title('Saccade delay frequency for no-stop trials','FontName','calibri','FontSize',12);
    set(gca,'Xlim',[0 600])
    set(gca,'TickDir','out','box','off');
    newpos =  get(gcf,'Position')/60;
    set(gcf,'PaperUnits','inches','PaperPosition',newpos);
    exportfigname=[cell2mat(regexp(directory,'\w+:\\\w+\\','match')),...
        'Analysis\Countermanding\',recname(1:end-4),'_PsyCurves'];
    %             print(psychoplots, '-dpng', '-noui', '-opengl','-r600', exportfigname);
    %            plot2svg([exportfigname,'.svg'],psychoplots, 'png');
    delete(psychoplots);
end

%% calculate SSRT
if ~(isempty(prevssds) || length(prevssds)==1)
    
    % test monotonicity and keep relevant inhibition function
    %     if (all(diff(probaresp)>=0) && length(probaresp)>=3) && ~all(diff(probaresp_diff)>=0)
    %         inhibfun=probaresp;
    %         ssds=narssdbins;
    %     else
    if all(diff(fullgauss_filtconv(probaresp_diff,1,1))>=0)
        inhibfun=probaresp_diff;
        if length(prevssds)>=2
            ssds=prevssds;
        else
            ssds=ssdvalues;
        end
        %     elseif (all(diff(probaresp)>=0) && length(probaresp)>=2) && ~all(diff(probaresp_diff)>=0)
        %         inhibfun=probaresp;
        %         ssds=narssdbins;
    else %monotonicity failure
        disp('failed to generate monotonic inhibition function')
        [overallMeanSSRT, inhibfun, ssds]=deal(NaN);
    end
    
    % calculate SSRT with Boucher et al's method
    try
        overallMeanSSRT= ...
            ssrt_bestfit(sacdelay.nsst', inhibfun', ssds', (mean(tachomc)+tachowidth));
    catch
        [overallMeanSSRT]=deal(NaN);
    end
else
    [overallMeanSSRT, inhibfun, ssds]=deal(NaN);
end

mssrt=overallMeanSSRT; %,meanIntSSRT,meanSSRT];
%     mssrt=round(nanmean(mssrt(mssrt>mean(tachomc)+10 & mssrt<130)));


if isnan(mssrt) || ~(mssrt>50 & mssrt<150) %get tachomc and lookup SSRT/tachomc fit. If fit missing, run SSRT_TachoMP
%     try
%         load([recname(1),'_tachoSSRTfit'],'fit');
%     catch
%         %SSRT_TachoMP
%     end

    % find reciprocal SSRT value
    try
        %get tacho curve midpoint
        tachomc=mean(tachomc);
        %estimate mssrt
%         mssrt=round(tachomc*fit.coeff(1)+fit.coeff(2));
%         if (tachomc>40 && tachomc<75)
             mssrt=round(tachomc+tachowidth);
%         elseif (tachomc>75 && tachomc<110)
%              mssrt=round(tachomc+tachowidth/2);
%         end
    catch

    end
end

% if mssrt < 50 ||  mssrt > 200
%    mssrt=NaN;
% end
%     
if ~(mssrt>70 & mssrt<130)
    userinfo=SetUserDir;
    load([userinfo.syncdir userinfo.slash 'behavData.mat']);
    SSRTs=[behavData.mssrt];
    if strcmp('R',recname(1))
        indivSSRTs=SSRTs(cellfun(@(x) strcmp(x,'Rigel'), {behavData.subject}));
    elseif strcmp('S',recname(1))
        indivSSRTs=SSRTs(cellfun(@(x) strcmp(x,'Sixx'), {behavData.subject}));
    elseif strcmp('H',recname(1))
        indivSSRTs=SSRTs(cellfun(@(x) strcmp(x,'Hilda'), {behavData.subject}));
    end
    indivSSRTs=indivSSRTs(~isnan(indivSSRTs));
    if ~isnan(mssrt)
        mssrt=round(mssrt/3+(mean(indivSSRTs)*2/3));
    else
        mssrt=mean(indivSSRTs);
    end
%     load([recname(1),'_evolSSRT'],'evolSSRT','foSSRT');
%     session=regexp(recname,'\d+','match');
%     if min(abs(evolSSRT(2,:)-str2num(session{1})))<=5
%         mssrt=round(mssrt/3+(evolSSRT(1,find(abs(evolSSRT(2,:)-str2num(session{1}))==min(abs(evolSSRT(2,:)-str2num(session{1}))),1)))*2/3);
%     else
%         mssrt=round(mssrt/3+foSSRT*2/3);
%     end
end

end
