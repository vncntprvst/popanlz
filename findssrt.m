function findssrt  %% calculate ssrt
%load('S99L2A5_13242.mat'); %old data format: random ssd
load('S116L4A6_15431.mat'); %new data format, four ssd : 100 150 200 250
%load('S117L4A6_12741.mat'); %new data format, four ssd : 100 150 200 250
%load('S118L4A5_13081.mat'); %new data format, four ssd : 80, 160, 210, 320.
load('R158L0A3_20310.mat'); % testing ssd 120,200,280,360, with Rigel: on first try he doesn't manage 280 and 360 
% for more analysis on delays and wrong trials, see gapstopdebug

%global allcodes alltimes;
% find stop trials
trialtypes=floor(allcodes(:,2)./10);
stoptrials=find(trialtypes==407);

stoptrialcodes=allcodes(stoptrials,:);
if find(stoptrialcodes(:,8)==1503,1) %recordings with benchmark code
    abortedstoptrials=~(floor(stoptrialcodes(:,9)./10)==507);
    noncancel=(stoptrialcodes(~abortedstoptrials,10)==16386) | ... % trials where a saccade is initiated 
    (stoptrialcodes(~abortedstoptrials,10)==17385);                % or subsequently breaking fixation
else
    abortedstoptrials=~(floor(stoptrialcodes(:,8)./10)==507);
    noncancel=(stoptrialcodes(~abortedstoptrials,9)==17385) | ...
    (stoptrialcodes(~abortedstoptrials,9)==16386);
end

% find delay times
stoptrialtimes=alltimes(stoptrials(~abortedstoptrials,:),:);
noncanceltimes=stoptrialtimes(noncancel,:);
canceltimes=stoptrialtimes(~noncancel,:);
noncancelssd=(noncanceltimes(:,8)-noncanceltimes(:,7))-3; % keep in mind that video synch adds ~55ms
cancelssd=(canceltimes(:,8)-canceltimes(:,7))-3;%3 ms added by the state transitions
ssdvalues=unique(cancelssd); 


% find withholding rate (ie, cancel) for each ssd
% allstoptrialcodes=stoptrialcodes(~abortedstoptrials,:);
% allnoncancel=allstoptrialcodes(noncancel,:);
% allcancel=allstoptrialcodes(~noncancel,:);

if length(ssdvalues)==4
withhldrate1=(length(find(cancelssd==ssdvalues(1)))./...
    (length(find(cancelssd==ssdvalues(1)))+length(find(noncancelssd==ssdvalues(1)))))*100
withhldrate2=(length(find(cancelssd==ssdvalues(2)))./...
    (length(find(cancelssd==ssdvalues(2)))+length(find(noncancelssd==ssdvalues(2)))))*100
withhldrate3=(length(find(cancelssd==ssdvalues(3)))./...
    (length(find(cancelssd==ssdvalues(3)))+length(find(noncancelssd==ssdvalues(3)))))*100
withhldrate4=(length(find(cancelssd==ssdvalues(4)))./...
    (length(find(cancelssd==ssdvalues(4)))+length(find(noncancelssd==ssdvalues(4)))))*100
else
    %collapse values
end



    % saccade delay for non-stop trials
    trialtot=length(datalign(1,1).allgreyareas);
    cuetimes=nan(trialtot,1);
    for trialnum=1:trialtot
    allgreytimes=datalign(1,1).allgreyareas{1,trialnum};
    cuetimes(trialnum)=allgreytimes(1,1);
    end
    sactime=datalign(1,1).alignidx;
    sacdelay=sactime-cuetimes;
    
    % calculating probabilities
    if size(datalign(end).ssd,2)==2
        %delaybincenters=unique(datalign(end).ssd(:,2))'+round(mean(datalign(end).ssd(:,1)-datalign(end).ssd(:,2)));
        [nccssdhist,delaybincenters]=hist(datalign(end).ssd(:,1),4);
        ccssdhist=hist(datalign(end-1).ssd(:,1),4);
    else
        delaybincenters=[163   213   273   313]; % Boucher et al values :[69   117   169   217]
        % adjust according to monkey history and actual distribution
        nccssdhist=hist(datalign(end).ssd(:,1),delaybincenters);
        ccssdhist=hist(datalign(end-1).ssd(:,1),delaybincenters);
    end

    probastop=nccssdhist./(nccssdhist+ccssdhist);
    
    %finding ssrt - first method
    
    %mean inhibition function
    mininhibfun=(sum((probastop(2:4)-probastop(1:3)).*delaybincenters(2:4)))./(max(probastop)-min(probastop));
    inhibfunssrt(1)=mean(sacdelay)-mininhibfun;
    % fitting weibull function
    if ~(length(find(probastop))<length(probastop))
        weibullcdfparam=wblfit(probastop);
        weibullcdf=wblcdf(probastop,weibullcdfparam(1),weibullcdfparam(2));
        %plot(weibullcdf)
        mininhibfunwb=(sum((weibullcdf(2:4)-weibullcdf(1:3)).*delaybincenters(2:4)))...
            ./(max(weibullcdf)-min(weibullcdf));
        inhibfunssrt(2)=mean(sacdelay)-mininhibfunwb;
    end
    
    %finding ssrt - second method
    delaydistribedges=min(sacdelay)-1:10:max(sacdelay)+1;
    delaydistribpdf=histc(sort(sacdelay),delaydistribedges);
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
    sufficientdelay=find(probastop(1)<=cumtrapz(delayfreq),1);
    intssrt(1)=delaydistribedges(sufficientdelay);
    sufficientdelay=find(probastop(2)<=cumtrapz(delayfreq),1);
    intssrt(2)=delaydistribedges(sufficientdelay);
    sufficientdelay=find(probastop(3)<=cumtrapz(delayfreq),1);
    intssrt(3)=delaydistribedges(sufficientdelay);
    sufficientdelay=find(probastop(4)<=cumtrapz(delayfreq),1);
    intssrt(4)=delaydistribedges(sufficientdelay);
    
    if inhibfunssrt>=(mean(intssrt)+std(intssrt))
        ssrt=mean(intssrt);
    else
        ssrt=mean([intssrt inhibfunssrt]);
    end
    
    % finally, adding SSRT to canceled trials alignment time, and selecting
    % latency matched no-stop trials
    
    alignnum=find(strcmp({datalign.alignlabel},'stop_cancel'));   
    ssd=datalign(1,alignnum).ssd(:,1);
    ssd1trials=ssd<ceil(delaybincenters(1));
    datalign(1,alignnum).alignidx=datalign(alignnum).alignidx+ssrt; %using a single value for ssrt
    datalign(1,alignnum).rasters=datalign(1,alignnum).rasters(ssd1trials,:);
    datalign(1,alignnum).trials=datalign(1,alignnum).trials(ssd1trials);
    datalign(1,alignnum).timefromtrig=datalign(1,alignnum).timefromtrig(ssd1trials);
    datalign(1,alignnum).timetotrig=datalign(1,alignnum).timetotrig(ssd1trials);
    datalign(1,alignnum).eyeh=datalign(1,alignnum).eyeh(ssd1trials,:);
    datalign(1,alignnum).eyev=datalign(1,alignnum).eyev(ssd1trials,:);
    datalign(1,alignnum).eyevel=datalign(1,alignnum).eyevel(ssd1trials,:);
    datalign(1,alignnum).allgreyareas=datalign(1,alignnum).allgreyareas(ssd1trials);
    datalign(1,alignnum).amplitudes=datalign(1,alignnum).amplitudes(ssd1trials);
    datalign(1,alignnum).peakvels=datalign(1,alignnum).peakvels(ssd1trials);
    datalign(1,alignnum).peakaccs=datalign(1,alignnum).peakaccs(ssd1trials);
    datalign(1,alignnum).bad=datalign(1,alignnum).bad(ssd1trials);
    
    %same restriction for non-cancel trials
    alignnum=find(strcmp({datalign.alignlabel},'stop_non_cancel'));   
    ssd=datalign(1,alignnum).ssd(:,1);
    ssd1trials=ssd<ceil(delaybincenters(1));
    datalign(1,alignnum).alignidx=datalign(alignnum).alignidx+ssrt; %using a single value for ssrt
    datalign(1,alignnum).rasters=datalign(1,alignnum).rasters(ssd1trials,:);
    datalign(1,alignnum).trials=datalign(1,alignnum).trials(ssd1trials);
    datalign(1,alignnum).timefromtrig=datalign(1,alignnum).timefromtrig(ssd1trials);
    datalign(1,alignnum).timetotrig=datalign(1,alignnum).timetotrig(ssd1trials);
    datalign(1,alignnum).eyeh=datalign(1,alignnum).eyeh(ssd1trials,:);
    datalign(1,alignnum).eyev=datalign(1,alignnum).eyev(ssd1trials,:);
    datalign(1,alignnum).eyevel=datalign(1,alignnum).eyevel(ssd1trials,:);
    datalign(1,alignnum).allgreyareas=datalign(1,alignnum).allgreyareas(ssd1trials);
    datalign(1,alignnum).amplitudes=datalign(1,alignnum).amplitudes(ssd1trials);
    datalign(1,alignnum).peakvels=datalign(1,alignnum).peakvels(ssd1trials);
    datalign(1,alignnum).peakaccs=datalign(1,alignnum).peakaccs(ssd1trials);
    datalign(1,alignnum).bad=datalign(1,alignnum).bad(ssd1trials);
    
    % last but not least, finding latency matched no-stop signal trials 
    
    alignnum=find(strcmp({datalign.alignlabel},'sac'));
    latmatch=sacdelay>=ceil(delaybincenters(1))+ssrt;
    
    datalign(1,alignnum).rasters=datalign(1,alignnum).rasters(latmatch,:);
    %          datalign(1,alignnum).alignidx=datalign(1,alignnum).alignidx; %    doesn't change
    datalign(1,alignnum).trials=datalign(1,alignnum).trials(latmatch);
    datalign(1,alignnum).timefromtrig=datalign(1,alignnum).timefromtrig(latmatch);
    datalign(1,alignnum).timetotrig=datalign(1,alignnum).timetotrig(latmatch);
    datalign(1,alignnum).eyeh=datalign(1,alignnum).eyeh(latmatch,:);
    datalign(1,alignnum).eyev=datalign(1,alignnum).eyev(latmatch,:);
    datalign(1,alignnum).eyevel=datalign(1,alignnum).eyevel(latmatch,:);
    datalign(1,alignnum).allgreyareas=datalign(1,alignnum).allgreyareas(latmatch);
    datalign(1,alignnum).amplitudes=datalign(1,alignnum).amplitudes(latmatch);
    datalign(1,alignnum).peakvels=datalign(1,alignnum).peakvels(latmatch);
    datalign(1,alignnum).peakaccs=datalign(1,alignnum).peakaccs(latmatch);
    datalign(1,alignnum).bad=datalign(1,alignnum).bad(latmatch);
    
end


