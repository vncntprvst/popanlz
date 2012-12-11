function findssrt  %% calculate ssrt

if strcmp(getenv('username'),'SommerVD')
    directory = 'C:\Data\Recordings\';
elseif  strcmp(getenv('username'),'DangerZone')
    directory = 'E:\Data\Recordings\';
else
    directory = 'B:\data\Recordings\';
end
slash = '\';

splitdataaligned=0;

%fname=('S99L2A5_13242.mat'); %old data format: random ssd
%fname=('S116L4A6_15431.mat'); %new data format, four ssd : 100 150 200 250
%fname=('S117L4A6_12741.mat'); %new data format, four ssd : 100 150 200 250
%fname=('S118L4A5_13081.mat'); %new data format, four ssd : 80, 160, 210, 320.
%fname=('R158L0A3_20310.mat'); % testing ssd 120,200,280,360, with Rigel: on first try he doesn't manage 280 and 360

%filestoload={'S99L2A5_13242.mat','S116L4A6_15431.mat','S117L4A6_12741.mat','S118L4A5_13081.mat'};

%% get gapstop filenames from procdata excel file
monknum=2;
    %Get number rows
    exl = actxserver('excel.application');
    exlWkbk = exl.Workbooks;
    exlFile = exlWkbk.Open([directory 'procdata.xlsx']);
    exlSheet = exlFile.Sheets.Item(monknum);% e.g.: 2 = Sixx
    robj = exlSheet.Columns.End(4);
    numrows = robj.row;
    % if numrows==1048576 %empty document
    %     numrows=1;
    % end
    Quit(exl);
    
    % read A (File Name) and G (Task) columns 
    [~,allfilenames] = xlsread([directory 'procdata.xlsx'],monknum,['A2:A' num2str(numrows)]);
    [~,alltasks] = xlsread([directory 'procdata.xlsx'],monknum,['G2:G' num2str(numrows)]);
    filestoload=allfilenames(strcmp(alltasks,'gapstop'));

%% pre-alloc
allnccssd=cell(length(filestoload),1);
allccssd=cell(length(filestoload),1);
allsacdelay=cell(length(filestoload),1);
alldlbincnt=cell(length(filestoload),1);

%% get data
for numfile=1:length(filestoload)
fname=filestoload{numfile};
try
load(fname);
load([fname,'_sac']); %dataaligned file
catch
    continue;
end
% for more analysis on delays and wrong trials, see gapstopdebug

%global allcodes alltimes;

%% find stop trials
trialtypes=floor(allcodes(:,2)./10);
stoptrials=find(trialtypes==407);

stoptrialcodes=allcodes(stoptrials,:);
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

%% find "desired" delay times (keep in mind that video synch adds ~55ms, so real ssd are variable)
% and calculate ssd if field doesn't exist
stoptrialtimes=alltimes(stoptrials(~abortedstoptrials,:),:);
noncanceltimes=stoptrialtimes(noncancel,:);
canceltimes=stoptrialtimes(~noncancel,:);

noncancelssd=(noncanceltimes(:,8)-noncanceltimes(:,7))-3; %3 ms added by the state transitions
cancelssd=(canceltimes(:,8)-canceltimes(:,7))-3; %3 ms added by the state transitions
if ~isfield(dataaligned,'ssd')
    if bmssd %benchmark
        nccssd=(noncanceltimes(:,9)-noncanceltimes(:,7))-3; %3 ms added by the state transitions
        ccssd=(canceltimes(:,9)-canceltimes(:,7))-3;
    else
        nccssd=noncancelssd;
        ccssd=cancelssd;
    end
end

ssdvalues=unique(cancelssd);


%% find withholding rate (ie, cancel) for each ssd
% allstoptrialcodes=stoptrialcodes(~abortedstoptrials,:);
% allnoncancel=allstoptrialcodes(noncancel,:);
% allcancel=allstoptrialcodes(~noncancel,:);

if length(ssdvalues)==4
    withhldrate1=(length(find(cancelssd==ssdvalues(1)))./...
        (length(find(cancelssd==ssdvalues(1)))+length(find(noncancelssd==ssdvalues(1)))))*100;
    withhldrate2=(length(find(cancelssd==ssdvalues(2)))./...
        (length(find(cancelssd==ssdvalues(2)))+length(find(noncancelssd==ssdvalues(2)))))*100;
    withhldrate3=(length(find(cancelssd==ssdvalues(3)))./...
        (length(find(cancelssd==ssdvalues(3)))+length(find(noncancelssd==ssdvalues(3)))))*100;
    withhldrate4=(length(find(cancelssd==ssdvalues(4)))./...
        (length(find(cancelssd==ssdvalues(4)))+length(find(noncancelssd==ssdvalues(4)))))*100;
else
    % forget it
end


%% saccade delay for non-stop trials

trialtot=length(dataaligned(1,1).allgreyareas);
cuetimes=nan(trialtot,1);
for trialnum=1:trialtot
    allgreytimes=dataaligned(1,1).allgreyareas{1,trialnum};
    cuetimes(trialnum)=allgreytimes(1,1);
end
sactime=dataaligned(1,1).alignidx;
sacdelay=sactime-cuetimes;

%% collecting data
if ~isfield(dataaligned,'ssd')
        allnccssd{numfile}=nccssd;
        allccssd{numfile}=ccssd;
        allsacdelay{numfile}=sacdelay;
        [~,delaybincenters]=hist(ccssd,4);
        alldlbincnt{numfile}=round(delaybincenters');
else
        allnccssd{numfile}=dataaligned(end).ssd(:,1);
        allccssd{numfile}=dataaligned(end-1).ssd(:,1);
    if size(dataaligned(end).ssd,2)==2
        [~,delaybincenters]=hist(dataaligned(end-1).ssd(:,1),4);
        alldlbincnt{numfile}=round(delaybincenters');
    else
        [~,bincnt]=hist(cancelssd,4);
        alldlbincnt{numfile}=round(bincnt');
    end
end

clearvars -except allnccssd allccssd allsacdelay alldlbincnt filestoload numfile splitdataaligned;
end

%% old code for single file processing
% if ~isfield(dataaligned,'ssd')
%         nccssdhist=hist(nccssd,4);
%         ccssdhist=hist(ccssd,4);
%         [~,delaybincenters]=hist(ccssd,4);
%         delaybincenters=round(delaybincenters);
% else
%     if size(dataaligned(end).ssd,2)==2
%         %delaybincenters=unique(dataaligned(end).ssd(:,2))'+round(mean(dataaligned(end).ssd(:,1)-dataaligned(end).ssd(:,2)));
%         nccssdhist=hist(dataaligned(end).ssd(:,1),4);
%         ccssdhist=hist(dataaligned(end-1).ssd(:,1),4);
%         [~,delaybincenters]=hist(dataaligned(end-1).ssd(:,1),4);
%     else
%         [~,bincnt]=hist(cancelssd,4);
%         delaybincenters=round(bincnt);
%         %delaybincenters=[163   213   273   313]; % Boucher et al values :[69   117   169   217]
%         % adjust according to monkey history and actual distribution
%         nccssdhist=hist(dataaligned(end).ssd(:,1),delaybincenters);
%         ccssdhist=hist(dataaligned(end-1).ssd(:,1),delaybincenters);
%     end
% end

%% calculating probabilities
global probaresp;
global narssdbins;

allnccssd=vertcat(allnccssd{:});
allccssd=vertcat(allccssd{:});
allsacdelay=vertcat(allsacdelay{:});
    %remove failed binings (cells that show [1,2,3,4])
    alldlbincnt=alldlbincnt(~cellfun(@(x) size(x,2)>size(x,1), alldlbincnt));
alldlbincnt=vertcat(alldlbincnt{:});

% sort unique delay bin centers
ssdbins=unique(sort(alldlbincnt));
% narrow ssds to those found more than once
narssdbins=ssdbins(hist(allnccssd,ssdbins)>1);
% bin canceled and non-canceled trials according to these ssds
nccssdhist=hist(allnccssd,narssdbins);
ccssdhist=hist(allccssd,narssdbins);
% find probability to respond
probaresp=nccssdhist'./(nccssdhist'+ccssdhist');
% figure
% plot(narssdbins,probaresp,'o','MarkerEdgeColor','k','MarkerFaceColor','c','MarkerSize',8)
% set(gca,'XTick',[1:24],'XTickLabel',narssdbins)
% Generalized linear model regression
regcoeffs = glmfit(narssdbins,[nccssdhist' (nccssdhist'+ccssdhist')],'binomial','link','probit');
respregfit = glmval(regcoeffs, narssdbins,'probit','size', (nccssdhist'+ccssdhist'));
figure
subplot(2,2,3)
plot(narssdbins, nccssdhist'./(nccssdhist'+ccssdhist'),'o',...
    narssdbins,respregfit./(nccssdhist'+ccssdhist'),'-','LineWidth',2)
xlabel('stop signal delay (ms)');
title('Probability to respond')

%% finding ssrt - first method: integration (constant SSRT)
delaydistribedges=min(allsacdelay)-1:10:max(allsacdelay)+1;
delaydistribhhist=histc(sort(allsacdelay),delaydistribedges);
emptyvalues=find(~delaydistribhhist);
if emptyvalues(1)==1
    emptyvalues=emptyvalues(2:end);
end
if emptyvalues(end)==length(delaydistribhhist)
    emptyvalues=emptyvalues(1:end-1);
end
for empt=1:length(emptyvalues)
    delaydistribhhist(emptyvalues(empt))=mean([delaydistribhhist(emptyvalues(empt)-1)...
        delaydistribhhist(emptyvalues(empt)+1)]);
end
%interpdelaydistribpdf=interp1(delaydistribedges,delaydistribpdf,min(totsadelay)-1:1:max(totsadelay)+1);
delayfreq=delaydistribhhist./sum(delaydistribhhist);
subplot(2,2,1)
plot(delaydistribedges,delayfreq,'-','LineWidth',2);
xlabel('Saccade delay (ms)');
title('Saccade delay frequency for no-stop trials')
subplot(2,2,2)
plot(delaydistribedges,cumtrapz(delayfreq),'-','LineWidth',2);
xlabel('Saccade delay (ms)');
title('Cumulative distribution of saccade delay')

%interpdelfrq=interpdelaydistribpdf./sum(interpdelaydistribpdf);
%     wbfitparam=wblfit(delayfreq);
%     delaywbfit=wblcdf(delayfreq,wbfitparam(1),wbfitparam(2));

% finding SSRT for each SSD
intssrt=nan(length(narssdbins),1);
for prbr=1:length(narssdbins)
sufficientdelay=find(probaresp(prbr)<=cumtrapz(delayfreq),1);
stopprocfinishline=delaydistribedges(sufficientdelay);
intssrt(prbr)=stopprocfinishline-narssdbins(prbr);
end


%% finding ssrt - second method (random SSRT)

%mean inhibition function
mininhibfun=(sum((probaresp(2:end)-probaresp(1:end-1)).*narssdbins(2:end)))./(max(probaresp)-min(probaresp));
inhibfunssrt(1)=mean(allsacdelay)-mininhibfun;
% fitting weibull function
if ~(length(find(probaresp))<length(probaresp))
    weibullcdfparam=wblfit(probaresp);
    weibullcdf=wblcdf(probaresp,weibullcdfparam(1),weibullcdfparam(2));
    %plot(weibullcdf)
    mininhibfunwb=(sum((weibullcdf(2:end)-weibullcdf(1:end-1)).*narssdbins(2:end)))...
        ./(max(weibullcdf)-min(weibullcdf));
    inhibfunssrt(2)=mean(allsacdelay)-mininhibfunwb;
end

%% averaging different estimates
if mean(inhibfunssrt)>=(mean(intssrt)+std(intssrt)) || mean(inhibfunssrt)<=(mean(intssrt)-std(intssrt))
    ssrt=mean(intssrt);
else
    %ssrt=mean([intssrt' inhibfunssrt]);
    ssrt=mean([mean(intssrt) inhibfunssrt])
end

subplot(2,2,4)
plot(narssdbins,intssrt,'-','LineWidth',2)
hold on
plot(narssdbins,ones(length(narssdbins),1)*mean(inhibfunssrt),'-g','LineWidth',2)
plot(narssdbins,ones(length(narssdbins),1)*ssrt,'-r','LineWidth',2)
legend({'constant SSRT method','random SSRT method','average different estimates'})
xlabel('stop signal delay (ms)');
title('Stop signal reaction time');

%% finally, adding SSRT to canceled trials alignment time, and selecting
% latency matched non-stop trials
if splitdataaligned
    alignnum=find(strcmp({dataaligned.alignlabel},'stop_cancel'));
    ssd=dataaligned(1,alignnum).ssd(:,1);
    ssd1trials=ssd<ceil(delaybincenters(1));
    dataaligned(1,alignnum).alignidx=dataaligned(alignnum).alignidx+ssrt; %using a single value for ssrt
    dataaligned(1,alignnum).rasters=dataaligned(1,alignnum).rasters(ssd1trials,:);
    dataaligned(1,alignnum).trials=dataaligned(1,alignnum).trials(ssd1trials);
    dataaligned(1,alignnum).timefromtrig=dataaligned(1,alignnum).timefromtrig(ssd1trials);
    dataaligned(1,alignnum).timetotrig=dataaligned(1,alignnum).timetotrig(ssd1trials);
    dataaligned(1,alignnum).eyeh=dataaligned(1,alignnum).eyeh(ssd1trials,:);
    dataaligned(1,alignnum).eyev=dataaligned(1,alignnum).eyev(ssd1trials,:);
    dataaligned(1,alignnum).eyevel=dataaligned(1,alignnum).eyevel(ssd1trials,:);
    dataaligned(1,alignnum).allgreyareas=dataaligned(1,alignnum).allgreyareas(ssd1trials);
    dataaligned(1,alignnum).amplitudes=dataaligned(1,alignnum).amplitudes(ssd1trials);
    dataaligned(1,alignnum).peakvels=dataaligned(1,alignnum).peakvels(ssd1trials);
    dataaligned(1,alignnum).peakaccs=dataaligned(1,alignnum).peakaccs(ssd1trials);
    dataaligned(1,alignnum).bad=dataaligned(1,alignnum).bad(ssd1trials);
    
    %same restriction for non-cancel trials
    alignnum=find(strcmp({dataaligned.alignlabel},'stop_non_cancel'));
    ssd=dataaligned(1,alignnum).ssd(:,1);
    ssd1trials=ssd<ceil(delaybincenters(1));
    dataaligned(1,alignnum).alignidx=dataaligned(alignnum).alignidx+ssrt; %using a single value for ssrt
    dataaligned(1,alignnum).rasters=dataaligned(1,alignnum).rasters(ssd1trials,:);
    dataaligned(1,alignnum).trials=dataaligned(1,alignnum).trials(ssd1trials);
    dataaligned(1,alignnum).timefromtrig=dataaligned(1,alignnum).timefromtrig(ssd1trials);
    dataaligned(1,alignnum).timetotrig=dataaligned(1,alignnum).timetotrig(ssd1trials);
    dataaligned(1,alignnum).eyeh=dataaligned(1,alignnum).eyeh(ssd1trials,:);
    dataaligned(1,alignnum).eyev=dataaligned(1,alignnum).eyev(ssd1trials,:);
    dataaligned(1,alignnum).eyevel=dataaligned(1,alignnum).eyevel(ssd1trials,:);
    dataaligned(1,alignnum).allgreyareas=dataaligned(1,alignnum).allgreyareas(ssd1trials);
    dataaligned(1,alignnum).amplitudes=dataaligned(1,alignnum).amplitudes(ssd1trials);
    dataaligned(1,alignnum).peakvels=dataaligned(1,alignnum).peakvels(ssd1trials);
    dataaligned(1,alignnum).peakaccs=dataaligned(1,alignnum).peakaccs(ssd1trials);
    dataaligned(1,alignnum).bad=dataaligned(1,alignnum).bad(ssd1trials);
    
    % last but not least, finding latency matched no-stop signal trials
    
    alignnum=find(strcmp({dataaligned.alignlabel},'sac'));
    latmatch=sacdelay>=ceil(delaybincenters(1))+ssrt;
    
    dataaligned(1,alignnum).rasters=dataaligned(1,alignnum).rasters(latmatch,:);
    %          dataaligned(1,alignnum).alignidx=dataaligned(1,alignnum).alignidx; %    doesn't change
    dataaligned(1,alignnum).trials=dataaligned(1,alignnum).trials(latmatch);
    dataaligned(1,alignnum).timefromtrig=dataaligned(1,alignnum).timefromtrig(latmatch);
    dataaligned(1,alignnum).timetotrig=dataaligned(1,alignnum).timetotrig(latmatch);
    dataaligned(1,alignnum).eyeh=dataaligned(1,alignnum).eyeh(latmatch,:);
    dataaligned(1,alignnum).eyev=dataaligned(1,alignnum).eyev(latmatch,:);
    dataaligned(1,alignnum).eyevel=dataaligned(1,alignnum).eyevel(latmatch,:);
    dataaligned(1,alignnum).allgreyareas=dataaligned(1,alignnum).allgreyareas(latmatch);
    dataaligned(1,alignnum).amplitudes=dataaligned(1,alignnum).amplitudes(latmatch);
    dataaligned(1,alignnum).peakvels=dataaligned(1,alignnum).peakvels(latmatch);
    dataaligned(1,alignnum).peakaccs=dataaligned(1,alignnum).peakaccs(latmatch);
    dataaligned(1,alignnum).bad=dataaligned(1,alignnum).bad(latmatch);
end
end


