recname='S145L2A4_9065';
subject={'Rigel','Sixx','Hilda'};

load(recname,'allcodes','alltimes','allbad');
load([recname,'_sac'],'dataaligned');

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

% calculate ssd if field doesn't exist
stoptrialtimes=alltimes(stoptrials(~abortedstoptrials,:),:);
noncanceltimes=stoptrialtimes(noncancel,:);
canceltimes=stoptrialtimes(~noncancel,:);

    if bmssd %benchmark
        nccssd=(noncanceltimes(:,9)-noncanceltimes(:,7))-3; %3 ms added by the state transitions
        ccssd=(canceltimes(:,9)-canceltimes(:,7))-3;
    else
        nccssd=(noncanceltimes(:,8)-noncanceltimes(:,7))-2; %2 ms added by the state transitions
        ccssd=(canceltimes(:,8)-canceltimes(:,7))-2;
    end

ssdvalues=unique([ccssd;nccssd]);

%% get saccade delay for signal-respond (noncancel) trial

    if bmssd %benchmark
        alldata.SRsacdelay=(noncanceltimes(:,10)-noncanceltimes(:,7))-6; %6 ms added by the state transitions
    else
        alldata.SRsacdelay=(noncanceltimes(:,9)-noncanceltimes(:,7))-6;
    end

%% saccade delay for non-stop trials
sacdataalign=find(arrayfun(@(x) strcmp(x.alignlabel,'sac'),dataaligned));
sacdelay=struct('all',[],'left',[],'right',[]);
for sacda=1:length(sacdataalign)
cuetimes=cellfun(@(x) x(1,1),dataaligned(1,sacdataalign(sacda)).allgreyareas);
sactime=dataaligned(1,sacdataalign(sacda)).alignidx;
    %left / right assymetry
    if sum(dataaligned(1,sacdataalign(sacda)).amplitudes)>0
        sacdelay.right=sactime-cuetimes; %anywhere to the right
    else
        sacdelay.left=sactime-cuetimes;  %anywhere to the left
    end
end
sacdelay.all=[sacdelay.left sacdelay.right];
alldata.NSSsacdelay=sacdelay;%[sacdelay{:}];

% NSS success rate
% NSS_targ=allcodes(floor(allcodes(:,7)./10)==684,8);
% alldata.NSSsuccessrate=sum(floor(NSS_targ./10)==704)/length(NSS_targ);

%% calculating probabilities for this session
[~,delaybincenters]=hist([ccssd;nccssd],4);
alldlbincnt=round(delaybincenters');
% sort unique delay bin centers
ssdbins=unique(sort(alldlbincnt));
% narrow ssds to those found more than once
    % narssdbins=ssdbins(hist(nccssd,ssdbins)>1); %not for
    % individual SSRT calculation
    narssdbins=ssdbins(hist(nccssd,ssdbins)>0);
    alldata.ssd=narssdbins;
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
alldata.inhibfun=probaresp;

if ~(isempty(narssdbins) || length(narssdbins)==1)

%% calculate SSRT with Boucher et al's method    
try
    [meanIntSSRT, meanSSRT, overallMeanSSRT]= ...
    ssrt_bestfit(sacdelay.all', probaresp', narssdbins');
catch
   sacdelay;
end
allmeanssrt=overallMeanSSRT;
alldata.meanIntSSRT=meanIntSSRT;
alldata.meanSSRT=meanSSRT;
alldata.overallMeanSSRT=overallMeanSSRT;
end
% clearvars -except nccssd ccssd allsacdelay alldlbincnt filestoload ...
%     numfile splitdataaligned allmeanssrt alldata andir emdirections subject ...
%     monknum colecalldata;

%% prepare plot
    CMdatfig=figure;
%     CMdatfigpos=get(CMdatfig,'Position');
    CMdatfigpos=[500 300 800 500];
    set(CMdatfig,'Position',CMdatfigpos,'Color','w');
    
        % subplot to display values
    hvaldisp = axes('Position', [.6, .6, .3, .3]);
    set(hvaldisp,'Visible','off');
        htitle=title('SSD ranges and SSRT values','FontName','calibri','FontSize',14);
        set(htitle,'Visible','on')

% fit sigmoid through inhibition function
        fitresult = sigmoidfit(narssdbins, probaresp);
        yfitval=fitresult(50:10:400); % This gives us the templates for the SSD range-dependant inhibition functions (six for Sixx, ha ha)

% plot overall inhibition function
subplot(1,2,1);
plot(narssdbins,probaresp,'Color',[0.25 0.25 0.25],'LineWidth',1.8);
hold on
plot([50:10:400],yfitval,'Color',[0.2 0.4 0.6],'LineStyle','-.','LineWidth',1.5);
title('Inhibition function','FontName','calibri','FontSize',15);
hxlabel=xlabel(gca,'Stop Signal Delays','FontName','calibri','FontSize',12);
set(gca,'Xlim',[50 400],'XTick',[100:50:350],'TickDir','out','box','off'); %'XTickLabel',[50:50:400]
hylabel=ylabel(gca,'P(Signal-respond)','FontName','calibri','FontSize',12);
set(gca,'Ylim',[0 1],'TickDir','out','box','off');

% print SSRT values
[meanIntSSRT, meanSSRT, overallMeanSSRT]= ...
    ssrt_bestfit(sacdelay.all', probaresp', narssdbins');   
text(0.5,0.6,['meanIntSSRT = ' num2str(round(meanIntSSRT))],'FontName','calibri','FontSize',10,'Parent',hvaldisp);
text(0.5,0.5,['meanSSRT = ' num2str(round(meanSSRT))],'FontName','calibri','FontSize',10,'Parent',hvaldisp);
text(0.5,0.4,['overallMeanSSRT = ' num2str(round(overallMeanSSRT))],'FontName','calibri','FontSize',10,'Parent',hvaldisp);
    % now with inhib function sigmoid fit
    [meanIntSSRT, meanSSRT, overallMeanSSRT]= ...
    ssrt_bestfit(sacdelay.all', yfitval', [50:10:400]); 
    text(0.5,0.3,['meanIntSSRT_s = ' num2str(round(meanIntSSRT))],'FontName','calibri','FontSize',10,'Parent',hvaldisp);
    text(0.5,0.2,['meanSSRT_s = ' num2str(round(meanSSRT))],'FontName','calibri','FontSize',10,'Parent',hvaldisp);
    text(0.5,0.1,['overallMeanSSRT_s = ' num2str(round(overallMeanSSRT))],'FontName','calibri','FontSize',10,'Parent',hvaldisp);


%prepare and plot saccade latency frequency
[saclatquant,saclatxlims]=hist(sacdelay.all,[0:25:500]);
saclatfreq=saclatquant./sum(saclatquant);

subplot(1,2,2);
plot(saclatxlims,saclatfreq,'Color','k','LineWidth',1.8);
title('Saccade Latency','FontName','calibri','FontSize',15);
hxlabel=xlabel(gca,'Saccade latency','FontName','calibri','FontSize',12);
set(gca,'Xlim',[0 500],'XTick',[0:50:500],'TickDir','out','box','off'); %'XTickLabel',[50:50:400]
hylabel=ylabel(gca,'Proportion','FontName','calibri','FontSize',12);
curylim=get(gca,'YLim');
set(gca,'Ylim',[0 curylim(2)],'TickDir','out','box','off'); 

%% export figure
exportfigname=['CMdat_',subject{2},'_',recname];
print(gcf, '-dpng', '-noui', '-opengl','-r600', exportfigname);
delete(CMdatfig);


