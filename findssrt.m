function findssrt  %% calculate ssrt

splitdataaligned=0;

%fname=('S99L2A5_13242.mat'); %old data format: random ssd
%fname=('S116L4A6_15431.mat'); %new data format, four ssd : 100 150 200 250
%fname=('S117L4A6_12741.mat'); %new data format, four ssd : 100 150 200 250
%fname=('S118L4A5_13081.mat'); %new data format, four ssd : 80, 160, 210, 320.
%fname=('R158L0A3_20310.mat'); % testing ssd 120,200,280,360, with Rigel: on first try he doesn't manage 280 and 360

filestoload={'S99L2A5_13242.mat','S116L4A6_15431.mat','S117L4A6_12741.mat','S118L4A5_13081.mat'};

allnccssd=[];
allccssd=[];
allsacdelay=[];
alldlbincnt=[];

for numfile=1:length(filestoload)
fname=filestoload{numfile};
load(fname);
load([fname(1:end-4),'_sac']); %dataaligned file
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
        allnccssd=[allnccssd; nccssd];
        allccssd=[allccssd; ccssd];
        allsacdelay=[allsacdelay;sacdelay];
        [~,delaybincenters]=hist(ccssd,4);
        alldlbincnt=[alldlbincnt;round(delaybincenters')];
else
        allnccssd=dataaligned(end).ssd(:,1);
        allccssd=dataaligned(end-1).ssd(:,1);
    if size(dataaligned(end).ssd,2)==2
        [~,delaybincenters]=hist(dataaligned(end-1).ssd(:,1),4);
        alldlbincnt=[alldlbincnt;round(delaybincenters')];
    else
        [~,bincnt]=hist(cancelssd,4);
        alldlbincnt=[alldlbincnt;round(bincnt')];
    end
end

clearvars -except allnccssd allccssd allsacdelay alldlbincnt filestoload numfile;
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
ssdbins=unique(sort(alldlbincnt));
narssdbins=ssdbins(hist(allnccssd,ssdbins)>1);
nccssdhist=hist(allnccssd,narssdbins);
ccssdhist=hist(allccssd,narssdbins);
probaresp=nccssdhist'./(nccssdhist'+ccssdhist');


%% finding ssrt - first method: integration (constant SSRT)
delaydistribedges=min(allsacdelay)-1:10:max(allsacdelay)+1;
delaydistribpdf=histc(sort(allsacdelay),delaydistribedges);
emptyvalues=find(~delaydistribpdf);
if emptyvalues(1)==1
    emptyvalues=emptyvalues(2:end);
end
if emptyvalues(end)==length(delaydistribpdf)
    emptyvalues=emptyvalues(1:end-1);
end
for empt=1:length(emptyvalues)
    delaydistribpdf(emptyvalues(empt))=mean([delaydistribpdf(emptyvalues(empt)-1)...
        delaydistribpdf(emptyvalues(empt)+1)]);
end
%interpdelaydistribpdf=interp1(delaydistribedges,delaydistribpdf,min(totsadelay)-1:1:max(totsadelay)+1);
delayfreq=delaydistribpdf./sum(delaydistribpdf);
%plot(cumtrapz(delayfreq));

%interpdelfrq=interpdelaydistribpdf./sum(interpdelaydistribpdf);
%     wbfitparam=wblfit(delayfreq);
%     delaywbfit=wblcdf(delayfreq,wbfitparam(1),wbfitparam(2));

% finding SSRT for each SSD
intssrt=nan(length(narssdbins),1);
for prbr=1:length(narssdbins)
sufficientdelay=find(probaresp(prbr)<=cumtrapz(delayfreq),1);
intssrt(prbr)=delaydistribedges(sufficientdelay);
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
    ssrt=mean([intssrt inhibfunssrt]);
end

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


