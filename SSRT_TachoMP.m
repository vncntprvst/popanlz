function [colecalldata,allnccssd,allccssd]=SSRT_TachoMP(monknum,emdirections)  %% calculate ssrt individually for multiple files, plus overall
% adapt following code to add calculation of tachometric curve and make
% SSRT/curve midpoint correlation

global probaresp;
global narssdbins;

% subject={'Rigel','Sixx','Hilda'};
%emdirections={'upward','up_left','leftward','down_left','downward','down_right','rightward','up_right','all'};
% emdirections={'all'};

if strcmp(getenv('username'),'SommerVD') || strcmp(getenv('username'),'vp35')
    directory = 'C:\Data\Recordings\';
elseif  strcmp(getenv('username'),'DangerZone')
    directory = 'E:\Data\Recordings\';
else
    directory = 'B:\data\Recordings\';
end
slash = '\';

cd(directory)
% cd([directory 'figures' slash 'ssrt']);

splitdataaligned=0;

%fname=('S99L2A5_13242.mat'); %old data format: random ssd
%fname=('S116L4A6_15431.mat'); %new data format, four ssd : 100 150 200 250
%fname=('S117L4A6_12741.mat'); %new data format, four ssd : 100 150 200 250
%fname=('S118L4A5_13081.mat'); %new data format, four ssd : 80, 160, 210, 320.
%fname=('R158L0A3_20310.mat'); % testing ssd 120,200,280,360, with Rigel: on first try he doesn't manage 280 and 360

%filestoload={'S99L2A5_13242.mat','S116L4A6_15431.mat','S117L4A6_12741.mat','S118L4A5_13081.mat'};

%% get gapstop filenames from procdata excel file
%monknum=2;

if monknum==1
    subject='Rigel';
elseif monknum==2
    subject='Sixx';
elseif monknum==3
    subject='Hilda';
end

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

%% run analysis for each direction then all together
for andir=1:length(emdirections)
    
    %% pre-alloc
    allnccssd=cell(length(filestoload),1);
    allccssd=cell(length(filestoload),1);
    allsacdelay=cell(length(filestoload),1);
    alldlbincnt=cell(length(filestoload),1);
    allmeanssrt=cell(length(filestoload),1);
    alldata=struct('dir',[],'NSSsacdelay',[],'NSSsuccessrate',[],'SRsacdelay',[],'SSDs',[],'meanRT',[],...
        'inhibfun',[],'t_inhibfun',[],'NRMSE',[],'monotonic',[],'meanIntSSRT',[],'meanSSRT',[],'overallMeanSSRT',[],...
        'intssrt',[],'simpleSSRT',[],'substrSSRT',[],'synthSSRT',[],'tacho',[],'xctr',[]);
    
 
    %% get data
    for numfile=1:length(filestoload)
        fname=filestoload{numfile};
        %select which file to use
        if size([which([eval('fname') '_REX.mat']);which([eval('fname') '_Sp2.mat'])],1)>1
            fname=[fname '_Sp2'];
        else
            fname=[fname '_REX'];
        end
        
        try
            load(fname,'allcodes','alltimes','allbad','saccadeInfo');
            try
                load([fname,'_sactgtrew'],'dataaligned');
                dataaligned=dataaligned{1}; % keep only sac alignment
            catch
                try
                    load([fname(1:end-4),'_sac'],'dataaligned');
                catch
                    continue
                end
            end
        catch
            continue;
        end
        % for more analysis on delays and wrong trials, see gapstopdebug
        
        %global allcodes alltimes;
        alldata(numfile).dir=emdirections(andir);
        
        %% get directions from eye movements traces
        % use code written for SummaryPlot to find directions
        for danum=1:length(dataaligned)
            sacdeg=nan(size(dataaligned(1,danum).trials,2),1);
            for eyetr=1:size(dataaligned(1,danum).trials,2)
                thissach=dataaligned(1,danum).eyeh(eyetr,dataaligned(1,danum).alignidx:dataaligned(1,danum).alignidx+100);
                thissacv=dataaligned(1,danum).eyev(eyetr,dataaligned(1,danum).alignidx:dataaligned(1,danum).alignidx+100);
                minwidth=5;
                [~, ~, thissacvel, ~, ~, ~] = cal_velacc(thissach,thissacv,minwidth);
                peakvel=find(thissacvel==max(thissacvel),1);
                sacendtime=peakvel+find(thissacvel(peakvel:end)<=...
                    (min(thissacvel(peakvel:end))+(max(thissacvel(peakvel:end))-min(thissacvel(peakvel:end)))/10),1);
                try
                    sacdeg(eyetr)=abs(atand((thissach(sacendtime)-thissach(1))/(thissacv(sacendtime)-thissacv(1))));
                catch
                    thissacv;
                end
                
                % sign adjustements
                if thissacv(sacendtime)<thissacv(1) % negative vertical amplitude -> vertical flip
                    sacdeg(eyetr)=180-sacdeg(eyetr);
                end
                if thissach(sacendtime)>thissach(1)%inverted signal: leftward is in postive range. Correcting to negative.
                    sacdeg(eyetr)=360-sacdeg(eyetr); % mirror image;
                end
            end
            % a quick fix to be able to put "upwards" directions together
            distrib=hist(sacdeg,3); %floor(length(sacdeg)/2)
            if max(bwlabel(distrib,4))>1 && distrib(1)>1 && distrib(end)>1 %=bimodal distribution with more than 1 outlier
                sacdeg=sacdeg+45;
                sacdeg(sacdeg>360)=-(360-(sacdeg(sacdeg>360)-45));
                sacdeg(sacdeg>0)= sacdeg(sacdeg>0)-45;
            end
            sacdeg=abs(median(sacdeg));
            
            if sacdeg>45/2 && sacdeg <= 45+45/2
                dataaligned(1,danum).dir='up_right';
            elseif sacdeg>45+45/2 && sacdeg <= 90+45/2
                dataaligned(1,danum).dir='rightward';
            elseif sacdeg>90+45/2 && sacdeg <= 135+45/2
                dataaligned(1,danum).dir='down_right';
            elseif sacdeg>135+45/2 && sacdeg < 180+45/2
                dataaligned(1,danum).dir='downward';
            elseif sacdeg>=180+45/2 && sacdeg <= 225+45/2
                dataaligned(1,danum).dir='down_left';
            elseif sacdeg>225+45/2 && sacdeg <= 270+45/2
                dataaligned(1,danum).dir='leftward';
            elseif sacdeg>270+45/2 && sacdeg <= 315+45/2
                dataaligned(1,danum).dir='up_left';
            else
                dataaligned(1,danum).dir='upward';
            end
        end
        if ~strcmp(emdirections(andir),'all')
            %check if file contains requested direction
            if ~(sum(cellfun(@(x) strcmp(x,emdirections(andir)), {dataaligned.dir})))
                %         dataaligned;
                continue;
            end
            % select directions
            %dataaligned=dataaligned(cellfun(@(x) (strcmp(x,'rightward') || strcmp(x,'leftward')), {dataaligned.dir}));
            dataaligned=dataaligned(cellfun(@(x) strcmp(x,emdirections(andir)), {dataaligned.dir}));
            
            %% find stop trials
            trialtypes=floor(allcodes(:,2));%./10 only /10 if pooling data together
            if length(emdirections)<9
                dirdig=find(cellfun(@(x) strcmp(x,emdirections(andir)), {'upward','up_left','leftward','down_left','downward','down_right','rightward','up_right','all'}))-1;
                stoptrials=find(trialtypes==(4070+dirdig));
            else
                stoptrials=find(trialtypes==(4070+andir-1)); %eg, 2 or 6 to 407 for left of right respectively
            end
            if isempty(stoptrials)
                %perhaps something went wrong in direction detection
                trialtypes;
                continue
            end
        else
            
            %% find stop trials
            trialtypes=floor(allcodes(:,2)./10);
            stoptrials=find(trialtypes==407);
        end
        
        %% get noncanceled trials
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
        
        %keep track of trial number for SSD transform
        goodstoptrials=stoptrials(~abortedstoptrials);
        noncanceltrials=goodstoptrials(noncancel);
        canceltrials=goodstoptrials(~noncancel);
        allsactrials=find(trialtypes==604);
        if find(stoptrialcodes(:,8)==1503,1)
            sactrials=allsactrials(floor(allcodes(allsactrials,9)./10)==704);
        else
            sactrials=allsactrials(floor(allcodes(allsactrials,8)./10)==704);
        end
        
        
        %% find delay times (keep in mind that video synch adds ~65ms, so real ssd are variable)
        
        % and calculate ssd if field doesn't exist
        stoptrialtimes=alltimes(stoptrials(~abortedstoptrials,:),:);
        noncanceltimes=stoptrialtimes(noncancel,:);
        canceltimes=stoptrialtimes(~noncancel,:);
        
        % if ~isfield(dataaligned,'ssd')
        if bmssd %benchmark
            nccssd=(noncanceltimes(:,9)-noncanceltimes(:,7))-3; %3 ms added by the state transitions
            ccssd=(canceltimes(:,9)-canceltimes(:,7))-3;
        else
            nccssd=(noncanceltimes(:,8)-noncanceltimes(:,7))-2; %2 ms added by the state transitions
            ccssd=(canceltimes(:,8)-canceltimes(:,7))-2;
        end
        
        % ssdvalues=unique(ccssd);
        
        %% get saccade delay for signal-respond (noncancel) trial
        
        if bmssd %benchmark
            alldata(numfile).SRsacdelay=(noncanceltimes(:,10)-noncanceltimes(:,7))-6; %6 ms added by the state transitions
        else
            alldata(numfile).SRsacdelay=(noncanceltimes(:,9)-noncanceltimes(:,7))-6;
        end
        
        %% find withholding rate (ie, cancel) for each ssd
        % allstoptrialcodes=stoptrialcodes(~abortedstoptrials,:);
        % allnoncancel=allstoptrialcodes(noncancel,:);
        % allcancel=allstoptrialcodes(~noncancel,:);
        
        % if length(ssdvalues)==4
        %     withhldrate1=(length(find(cancelssd==ssdvalues(1)))./...
        %         (length(find(cancelssd==ssdvalues(1)))+length(find(noncancelssd==ssdvalues(1)))))*100;
        %     withhldrate2=(length(find(cancelssd==ssdvalues(2)))./...
        %         (length(find(cancelssd==ssdvalues(2)))+length(find(noncancelssd==ssdvalues(2)))))*100;
        %     withhldrate3=(length(find(cancelssd==ssdvalues(3)))./...
        %         (length(find(cancelssd==ssdvalues(3)))+length(find(noncancelssd==ssdvalues(3)))))*100;
        %     withhldrate4=(length(find(cancelssd==ssdvalues(4)))./...
        %         (length(find(cancelssd==ssdvalues(4)))+length(find(noncancelssd==ssdvalues(4)))))*100;
        % else
        %     % forget it
        % end
        
        
        %% saccade delay for non-stop trials
        % first method: separate left and right
        sacdataalign=find(arrayfun(@(x) strcmp(x.alignlabel,'sac'),dataaligned));
        sacdelay=struct('all',[],'left',[],'right',[]);
        for sacda=1:length(sacdataalign)
            try
            cuetimes=cellfun(@(x) x(1,1),dataaligned(1,sacdataalign(sacda)).allgreyareas);
            catch
                continue
            end
            sactime=dataaligned(1,sacdataalign(sacda)).alignidx;
            %left / right assymetry
            if sum(dataaligned(1,sacdataalign(sacda)).amplitudes)>0
                sacdelay.right=sactime-cuetimes; %anywhere to the right
            else
                sacdelay.left=sactime-cuetimes;  %anywhere to the left
            end
        end
        % sacdelay.all=[sacdelay.left sacdelay.right];
        
        % second method: all good saccade from non-stop trials (may yield slightly
        % different results than with left/right parsing method above
        % check with: sactrials(~ismember(sactrials,sort([dataaligned(1,1).trials dataaligned(1,2).trials])))
        % and (after method below): sactrials(~ismember(sactrials,find(sum(allgoodsacs,2))))
        
        alllats=reshape({saccadeInfo.latency},size(saccadeInfo));
        alllats=alllats';%needs to be transposed because the logical indexing below will be done column by column, not row by row
        allgoodsacs=~cellfun('isempty',reshape({saccadeInfo.latency},size(saccadeInfo)));
        %weeding out bad trials that are not stop trials
        allgoodsacs(logical(allbad)'&trialtypes~=407,:)=0;
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
        
        alldata(numfile).NSSsacdelay=sacdelay;%[sacdelay{:}];
        
        % adjust sactrials
        sactrials=sactrials(ismember(sactrials,find(sum(allgoodsacs,2))));
        
        trialtot=length(dataaligned(1,1).allgreyareas);
        cuetimes=nan(trialtot,1);
        for trialnum=1:trialtot
            allgreytimes=dataaligned(1,1).allgreyareas{1,trialnum};
            cuetimes(trialnum)=allgreytimes(1,1);
        end
        sactime=dataaligned(1,1).alignidx;
        % sacdelay=sactime-cuetimes;
        
        % NSS success rate
        NSS_targ=allcodes(floor(allcodes(:,7)./10)==684,8);
        alldata(numfile).NSSsuccessrate=sum(floor(NSS_targ./10)==704)/length(NSS_targ);
        
        %% comparing probability density of RTs on trials with or without stop signal
        % alldata(numfile).SRsacdelay % stop signal trials
        % nccssd % corresponding ssd
        % sacdelay.all % no-stop signal trials
        [saclatquant,saclatxlims]=hist(sacdelay.all,[0:50:max([max(sacdelay.all) 500])]);
        saclatfreq=saclatquant./sum(saclatquant);
        %unique ssds index
        [~,~,ussdidx]=unique(nccssd);
        %allocate
        numdpperssd=nan(1,max(ussdidx));
        simpleSSRT=nan(1,max(ussdidx));
        % plots
        % figure;
        % cc=lines(max(ussdidx))+0.25;
        % cc(cc>1)=1;
        % plot(saclatxlims,saclatfreq,'Color','k','LineWidth',1.8);
        % legendstr{1}=('no-stop signal')
        % hold on
        for plotrt=1:max(ussdidx)
            srsaclat=alldata(numfile).SRsacdelay(ussdidx==plotrt);
            numdpperssd(plotrt)=length(srsaclat);
            srsaclatfreq=hist(srsaclat,length(srsaclat))./length(srsaclat);
            %     plot(srsaclat,srsaclatfreq,'Color',cc(plotrt,:),'LineWidth',1.8);
            %     legendstr{plotrt+1}=(['ssd= ' num2str(nccssd(find(ussdidx==plotrt,1)))]);
            
            % substract SSD from earlist saccade latencies of non-cenceled trials
            if length(srsaclat)>1
                simpleSSRT(plotrt)=round(mean(srsaclat)-nccssd(find(ussdidx==plotrt,1)));
            end
        end
        alldata(numfile).simpleSSRT=simpleSSRT;
        % legend(legendstr{:})
        % set(gca,'ylim',[0 1.1]);
        
        %% collecting data
        if ~isfield(dataaligned,'ssd')
            allnccssd{numfile}=nccssd;
            allccssd{numfile}=ccssd;
            allsacdelay{numfile}=sacdelay.all;
            %         if length(ccssd)>4
            %         adjust delay bin centers to number of data par ssd value
            [~,delaybincenters]=hist([ccssd;nccssd],round(median(numdpperssd)));
            alldlbincnt{numfile}=round(delaybincenters');
            %         else
            %             alldlbincnt{numfile}=ccssd;
            %         end
        else
            allnccssd{numfile}=dataaligned(end).ssd(:,1);
            allccssd{numfile}=dataaligned(end-1).ssd(:,1);
            %alldata(numfile).NSSsacdelay=dataaligned(1).sacdelay(:,1);
            if size(dataaligned(end).ssd,2)==2
                [~,delaybincenters]=hist(dataaligned(end-1).ssd(:,1),4);
                alldlbincnt{numfile}=round(delaybincenters');
            else
                [~,bincnt]=hist(ccssd,4);
                alldlbincnt{numfile}=round(bincnt');
            end
        end
        
        %% calculating probabilities for this session
        % sort unique delay bin centers
        ssdbins=unique(sort(alldlbincnt{numfile}));
        % narrow ssds to those found more than once
        % narssdbins=ssdbins(hist(allnccssd{numfile},ssdbins)>1);
        % not for individual SSRT calculation!!
        narssdbins=ssdbins;
        % bin canceled and non-canceled trials according to these ssds
        try
            nccssdhist=hist(allnccssd{numfile},narssdbins);
            ccssdhist=hist(allccssd{numfile},narssdbins);
            % find probability to respond
            probaresp=(floor((nccssdhist'./(nccssdhist'+ccssdhist')).*10))/10;
        catch
            probaresp=1;
            narssdbins=0;
        end
        alldata(numfile).inhibfun=probaresp;
        %% transform SSD: adjust relative to RT (see Nelson et al 2010)
        meanRT=mean(sacdelay.all);
        t_nccssd=nccssd;
        t_ccssd=ccssd;
        for transsd=1:length(noncanceltrials)
            %mean RT value around stop signal trial ("epoch RT")
                epochRT=nanmean([sacdelay.all(find(sactrials<noncanceltrials(transsd),1,'last')),...
                sacdelay.all(find(sactrials>noncanceltrials(transsd),1))]);
            %transform
            try
                t_nccssd(transsd)=t_nccssd(transsd)-round(epochRT-meanRT);
            catch
                continue
            end
        end
        for transsd=1:length(canceltrials)
            %mean RT value around stop signal trial ("epoch RT")
                epochRT=mean([sacdelay.all(find(sactrials<canceltrials(transsd),1,'last')),...
                sacdelay.all(find(sactrials>canceltrials(transsd),1))]);
            %transform
            try
            t_ccssd(transsd)=t_ccssd(transsd)-round(epochRT-meanRT);
            catch
                continue
            end
        end
        if length(t_ccssd)>4
            [~,t_delaybincenters]=hist([t_ccssd;t_nccssd],4);
            t_ssdbins=round(t_delaybincenters');
        else
            t_ssdbins=unique(sort(t_ccssd));
        end
        try
            t_nccssdhist=hist(t_nccssd,t_ssdbins);
            t_ccssdhist=hist(t_ccssd,t_ssdbins);
            t_probaresp=(floor((t_nccssdhist'./(t_nccssdhist'+t_ccssdhist')).*10))/10;
        catch
            t_probaresp=1;
            t_ssdbins=0;
        end
        alldata(numfile).t_inhibfun=t_probaresp;
        
        %% get tachometric curve. Build matrix with 1. SSDs 2. RTs 3. success
        SSDRTs=inf(length(allbad),3); % may differ from sum(~allbad)
        ccssd(ismember(ccssd,unique(ccssd-1)))=ccssd(ismember(ccssd,unique(ccssd-1)))+1;
        SSDRTs(canceltrials,1)=ccssd;
        nccssd(ismember(nccssd,unique(nccssd-1)))=nccssd(ismember(nccssd,unique(nccssd-1)))+1;
        SSDRTs(logical(sum(allncsacs,2)),1)=nccssd(ismember(goodstoptrials(noncancel),find(sum(allncsacs,2)))); %got to remove some ssd when sacs are not available
        SSDRTs(logical(sum(allgoodsacs,2)),2)=allsacdelay{numfile};
        SSDRTs(logical(sum(allncsacs,2)),2)=ncsacdelay;
        SSDRTs(logical(sum(allncsacs,2)),3)=0;
        SSDRTs(logical(allbad),3)=0;
        SSDRTs([find(sum(allgoodsacs,2));canceltrials],3)=1;
        try
            [xctr xtach tach rPTc rPTe] = tachCM2(SSDRTs);
        catch
            [xctr xtach tach rPTc rPTe] = deal(NaN);
        end
        alldata(numfile).tacho=tach;
        alldata(numfile).xctr=xctr;
        
        if ~(isempty(narssdbins) || length(narssdbins)==1)
            % % try to narrow inhibition function to monotonic part (if max value is
            % % reached before longest SSD)
            % if find(probaresp==(max(probaresp)))~=length(probaresp)
            %     if (max(probaresp)/probaresp(end))>1.1 %then it's really skewed
            %         narssdbins=narssdbins(1:find(probaresp==(max(probaresp))));
            %         probaresp=probaresp(1:find(probaresp==(max(probaresp))));
            %     end
            % end
            
            %% test monotonicity
            if all(diff(probaresp)>=0) || all(diff(t_probaresp)>=0)
                alldata(numfile).monotonic=1;
                
                %%plot inhibition function and sigmoid
                % plot mean of actual data
                %         IFploth=figure;
                %         plot(narssdbins,probaresp,'LineWidth',1.5);
                %         hold on
                %         plot(t_ssdbins,t_probaresp,'LineWidth',1.5,'Color','green');
                % calculate sigmoid fit
                %         fittime=0:10:500;
                %         [fitresult,gof] = sigmoidfit(narssdbins, probaresp);
                % the goodness of fit contains
                % sse Sum of squares due to error
                % R2  Coefficient of determination
                % adjustedR2  Degree-of-freedom adjusted coefficient of determination
                % stdError    Root mean squared error (standard error)
                
                %         yfitval=fitresult(fittime); % This gives us the templates for the SSD range-dependant inhibition functions (six for Sixx, ha ha)
                
                %         try
                %         t_fitresult = sigmoidfit(t_ssdbins, t_probaresp);
                %         t_yfitval=t_fitresult(fittime);
                %         catch
                %         end
                %
                
                
                %         % do some alignement between sigmoid and actual data
                %         ifstart=fittime(find(round(yfitval{bins}.*100)./100>=mean(fdsplitif{bins}(1,:),2),1)); %finds when the sig has same value has mean first inhibfun value
                %         ssdstart=fittime(find(fittime>mean(fdsplitssd{bins}(1,:),2),1));% fnd when mean first ssd value occurs
                %         sigdiff=ifstart-ssdstart;
                %         try
                %         fittime=fittime-sigdiff; %shift time
                %         catch
                %             sigdiff
                %         end
                % plot sigmoid
                %         plot(fittime(fittime>0),yfitval(fittime>0),'Color','blue','LineStyle','.','LineWidth',0.8);
                %         plot(fittime(fittime>0),t_yfitval(fittime>0),'Color','green','LineStyle','.','LineWidth',0.8);
                
                %% get residuals
                %From the stdError we can calculate the normalized root-mean-square deviation or error (NRMSD or NRMSE),...
                %that is the RMSD divided by the range of observed values of a variable being predicted.
                %         NRMSE=gof.rmse/(max(probaresp)-min(probaresp));
                %         t_NRMSE=t_gof.rmse/(max(t_probaresp)-min(t_probaresp));
                %         %The value is often expressed as a percentage, where lower values indicate less residual variance.
                %         alldata(numfile).NRMSE=round(NRMSE*100);
                
                %% print values
                
                %         text(300,0.9,['Residual variance is ',num2str(alldata(numfile).NRMSE),char(37)]);
                %         text(280,0.8,[num2str(t_NRMSE),char(37)]);
                %         text(300,0.7,['IF transform monotonic ',num2str(alldata(numfile).monotonic)]);
                %         text(280,0.6,num2str(all(diff(t_probaresp)>=0)));
                
                %% calculate SSRT with Boucher et al's method
                try
                    if all(diff(probaresp)>=0)
                        [meanIntSSRT, meanSSRT, overallMeanSSRT]= ...
                            ssrt_bestfit(sacdelay.all', probaresp', narssdbins');
                    else
                        [meanIntSSRT, meanSSRT, overallMeanSSRT]= ...
                            ssrt_bestfit(sacdelay.all', t_probaresp(t_ssdbins>0)', t_ssdbins(t_ssdbins>0)');
                    end
                    alldata(numfile).meanIntSSRT=meanIntSSRT;
                    alldata(numfile).meanSSRT=meanSSRT;
                    alldata(numfile).overallMeanSSRT=overallMeanSSRT;
                    
                    %% calculate SSRT with home-made code
                    % finding ssrt - first method: integration (constant SSRT)
                    delaydistribedges=min(sacdelay.all)-1:10:max(sacdelay.all)+1;
                    delaydistribhhist=histc(sort(sacdelay.all),delaydistribedges);
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
                    % subplot(2,2,1)
                    % plot(delaydistribedges,delayfreq,'-b','LineWidth',2);
                    % xlabel('Saccade delay (ms)');
                    % title('Saccade delay frequency for no-stop trials')
                    % subplot(2,2,2)
                    % plot(delaydistribedges,cumtrapz(delayfreq),'-','LineWidth',2);
                    % xlabel('Saccade delay (ms)');
                    % title('Cumulative distribution of saccade delay')
                    
                    %interpdelfrq=interpdelaydistribpdf./sum(interpdelaydistribpdf);
                    %     wbfitparam=wblfit(delayfreq);
                    %     delaywbfit=wblcdf(delayfreq,wbfitparam(1),wbfitparam(2));
                    
                    % finding SSRT for each SSD
                    if all(diff(probaresp)>=0)
                        intssrt=nan(length(narssdbins),1);
                        for prbr=1:length(narssdbins)
                            sufficientdelay=find(probaresp(prbr)<=cumtrapz(delayfreq),1);
                            stopprocfinishline=delaydistribedges(sufficientdelay);
                            intssrt(prbr)=stopprocfinishline-narssdbins(prbr);
                        end
                    else
                        t_probaresp=t_probaresp(t_ssdbins>0);
                        t_ssdbins=t_ssdbins(t_ssdbins>0);
                        intssrt=nan(length(t_ssdbins),1);
                        for prbr=1:length(t_ssdbins)
                            sufficientdelay=find(t_probaresp(prbr)<=cumtrapz(delayfreq),1);
                            stopprocfinishline=delaydistribedges(sufficientdelay);
                            intssrt(prbr)=stopprocfinishline-t_ssdbins(prbr);
                        end
                    end
                    alldata(numfile).intssrt=intssrt;
                    
                    %% simple substraction method: SSRT is the difference between the mean RT
                    % on no-signal trials and the midpoint of the inhibition function,
                    theta=atan((probaresp(2)-probaresp(1))/(narssdbins(2)-narssdbins(1)));
                    midpoint=round(cos(theta)*((0.5-probaresp(1))/sin(theta))+narssdbins(1));
                    alldata(numfile).substrSSRT=mean(sacdelay.all)-midpoint;
                    
                    % new synthetic SSRT
                        synthSSRT=[alldata(numfile).substrSSRT round(nanmean(simpleSSRT)) round(nanmean(intssrt))];
                        try
                            synthSSRT=mean(synthSSRT(synthSSRT>40 & synthSSRT<150));
                        catch
                            synthSSRT=NaN;
                        end
                        alldata(numfile).synthSSRT=synthSSRT;
                    
%                     if ~isempty(nanmean(simpleSSRT)) && (nanmean(simpleSSRT)>50 && nanmean(simpleSSRT)<150)
%                         xctr;
%                         tach;
%                         xtach;
%                     end
                    
                    
                    % keep values
                    alldata(numfile).meanRT=meanRT;
                    if all(diff(probaresp)>=0)
                        alldata(numfile).SSDs=narssdbins';
                    else
                        alldata(numfile).SSDs=t_ssdbins';
                    end
                    
                    if abs(diff([mean(intssrt),meanIntSSRT]))<10
                        
                        %% print more values
                        %         text(280,0.9,['meanRT ',num2str(meanRT)]);
                        %         text(280,0.8,['SSDs ',num2str(narssdbins')]);
                        %         text(280,0.7,['meanIntSSRT ',num2str(meanIntSSRT)]);
                        %         text(280,0.6,['intssrt ', num2str(intssrt')]);
                        %         text(280,0.5,['overall SSRT ', num2str(mean([intssrt; meanIntSSRT]))]);
                        %
                        %         hold off;
                        %         close(IFploth);
                        
%                         alldata(numfile).overallMeanSSRT=mean([mean(intssrt) meanIntSSRT]);
                        
                    end
                catch
                    sacdelay;
                end
                
            else
                alldata(numfile).monotonic=0;
            end
            
            %% finding ssrt - second method (random SSRT)
            
            % %mean inhibition function
            % mininhibfun=(sum((probaresp(2:end)-probaresp(1:end-1)).*narssdbins(2:end)))./(max(probaresp)-min(probaresp));
            % inhibfunssrt(1)=mean(sacdelay.all)-mininhibfun;
            % % fitting weibull function
            % if ~(length(find(probaresp))<length(probaresp))
            %     weibullcdfparam=wblfit(probaresp);
            %     weibullcdf=wblcdf(probaresp,weibullcdfparam(1),weibullcdfparam(2));
            %     %plot(weibullcdf)
            %     mininhibfunwb=(sum((weibullcdf(2:end)-weibullcdf(1:end-1)).*narssdbins(2:end)))...
            %         ./(max(weibullcdf)-min(weibullcdf));
            %     inhibfunssrt(2)=mean(sacdelay.all)-mininhibfunwb;
            % end
            %
            % %% averaging different estimates
            % if mean(inhibfunssrt)>=(mean(intssrt)+std(intssrt)) || mean(inhibfunssrt)<=(mean(intssrt)-std(intssrt))
            %     ssrt=mean(intssrt);
            % else
            %     %ssrt=mean([intssrt' inhibfunssrt]);
            %     ssrt=mean([mean(intssrt) inhibfunssrt]);
            % end
            
            % allmeanssrt{numfile}=overallMeanSSRT;
            % alldata(numfile).meanIntSSRT=meanIntSSRT;
            % alldata(numfile).meanSSRT=meanSSRT;
            % alldata(numfile).overallMeanSSRT=overallMeanSSRT;
            
        end
        clearvars -except allnccssd allccssd allsacdelay alldlbincnt filestoload ...
            numfile splitdataaligned allmeanssrt alldata andir emdirections subject ...
            monknum colecalldata;
    end
    
    %% plot evolution of SSRT with respect to SSDs and RT
SSRTs=cell(7,1);
    SSRTs{1,:}={alldata.substrSSRT};
    SSRTs{1,:}=cellfun(@(x) round(x), SSRTs{1,:},'UniformOutput',false);
    [SSRTs{1,:}{cellfun('isempty',SSRTs{1,:})}]=deal(0);
    SSRTs{2,:}={alldata.meanSSRT};
    SSRTs{2,:}=cellfun(@(x) round(x), SSRTs{2,:},'UniformOutput',false);
    [SSRTs{2,:}{cellfun('isempty',SSRTs{2,:})}]=deal(0);
    SSRTs{3,:}={alldata.meanIntSSRT};
    SSRTs{3,:}=cellfun(@(x) round(x), SSRTs{3,:},'UniformOutput',false);
    [SSRTs{3,:}{cellfun('isempty',SSRTs{3,:})}]=deal(0);
     SSRTs{4,:}={alldata.intssrt};
    SSRTs{4,:}=cellfun(@(x) round(nanmean(x)), SSRTs{4,:},'UniformOutput',false);
    [SSRTs{4,:}{cellfun('isempty',SSRTs{4,:})}]=deal(0);
    SSRTs{5,:}={alldata.simpleSSRT};
    SSRTs{5,:}=cellfun(@(x) round(nanmean(x)), SSRTs{5,:},'UniformOutput',false);
    [SSRTs{5,:}{cellfun('isempty',SSRTs{5,:})}]=deal(0);
    SSRTs{6,:}={alldata.overallMeanSSRT};
    SSRTs{6,:}=cellfun(@(x) round(x), SSRTs{6,:},'UniformOutput',false);
    [SSRTs{6,:}{cellfun('isempty',SSRTs{6,:})}]=deal(0);
    SSRTs{7,:}={alldata.synthSSRT};
    SSRTs{7,:}=cellfun(@(x) round(nanmean(x)), SSRTs{7,:},'UniformOutput',false);
    [SSRTs{7,:}{cellfun('isempty',SSRTs{7,:})}]=deal(0);
    %look at corresponding tacho midcurve value
    alltachomc={alldata.xctr};
    alltachomc=cellfun(@(x) round(nanmean(x)), alltachomc,'UniformOutput',false);
    [alltachomc{cellfun('isempty',alltachomc)}]=deal(0);

    %look at correlation coefficients
    figure(22)
    hold on
    figure(21);
    SSRTs_idx=cell(7,1);
    tachomc_idx=[alltachomc{:}]>20 & [alltachomc{:}]<120; 
    [tachomcdata,tachomcdata_idx]=sort([alltachomc{tachomc_idx}]);
    plot(tachomcdata,'LineWidth',1.5,'Color','k')
    hold on
%     RTtach={alldata(tachomc_idx).meanRT}; 
%     RTtach=RTtach(tachomcdata_idx)
%     plot([RTtach{~cellfun('isempty',RTtach)}],tachomcdata(~cellfun('isempty',RTtach)))
%     % no correlation detected between RT and tachomc
cc=lines(size(SSRTs,1));
    for ssrtidxnum=1:7
        SSRTs_idx{ssrtidxnum,:}=[SSRTs{ssrtidxnum,:}{:}]>40 & [SSRTs{ssrtidxnum,:}{:}]<150;
        [SSRTdata,~]=sort([SSRTs{ssrtidxnum,:}{SSRTs_idx{ssrtidxnum,:}}]);
        if ~isempty(SSRTdata)
        figure(21)
        plot(SSRTdata,'Color',cc(ssrtidxnum,:),'LineWidth',1.5)
        figure(22)
        tachosssrtidx=tachomc_idx & SSRTs_idx{ssrtidxnum,:};  
        scatter([alltachomc{tachosssrtidx}],[SSRTs{ssrtidxnum,:}{tachosssrtidx}],50,cc(ssrtidxnum,:))
        end
    end
    figure(21)
        legend('tachomcdata','substrSSRT','intssrt','simpleSSRT','synthSSRT','Location','SouthEast')
        
%% evolution of SSRT over time
    % find which SSRTs have appropriate values, and keep mean value
    nmsssrt=cell(size(SSRTs));
    for apssrt=1:size(SSRTs,1)
        if sum([SSRTs{apssrt,:}{:}]>40 & [SSRTs{apssrt,:}{:}]<150)>0
          nmsssrt{apssrt}=nanmean([SSRTs{apssrt,:}{[SSRTs{apssrt,:}{:}]>40 & [SSRTs{apssrt,:}{:}]<150}]);
        end
    end
    
    % final, overall value
    foSSRT=round(nanmean([nmsssrt{:}]));
    
    % evolution over time (row 1: SSRT, row 2: session)
    evolSSRT(1,:)=round(([SSRTs{7,:}{[SSRTs{7,:}{:}]>40 & [SSRTs{7,:}{:}]<150}])/2)+round(foSSRT/2);
    evolSSRT(2,:)=find([SSRTs{7,:}{:}]>40 & [SSRTs{7,:}{:}]<150);
    
moreplots=0
if moreplots

    plot(fidx(nsidx),meanvalsevol);
    set(gca,'ylim',[0 150]);
    xlabel('recording session')
    ylabel('SSRT (ms)');
    title('evolution of SSRTs over time');
    
    gsidx=(SSRTs>50 & SSRTs<150);
    SSRTs=SSRTs(gsidx); % round(mean(SSRTs)) Rigel=89;
    meanRTs=round([alldata(fidx(gsidx)).meanRT]);
    maxSSDs=cellfun(@(x) max(x), {alldata(fidx(gsidx)).SSDs});
    
    
    figure;
    plot(maxSSDs,'Color','green','LineWidth',1.5);
    legend('SSDs')
    hold on
    axes
    plot(meanRTs,'Color','red','LineWidth',1.5);
    legend('RTs')
    set(gca,'Color','none')
    axes
    plot(SSRTs,'Color','blue','LineWidth',1.5);
    set(gca,'Color','none')
    legend('SSRTs')
    
    figure;
    subplot(3,1,1)
    [smeanRTs,sidx]=sort(meanRTs);
    plot(smeanRTs,SSRTs(sidx));
    xlabel('mean RT (ms)')
    ylabel('SSRT (ms)')
    title('SSRT value as a function of RT')
    subplot(3,1,2)
    plot(fidx(gsidx),SSRTs);
    xlabel('recording session')
    ylabel('SSRT (ms)');
    title('evolution of SSRTs over time');
    subplot(3,1,3)
    plot(fidx(gsidx),SSRTs-meanRTs)
    title('evolution of SSRT - mean RT')
    xlabel('recording session')
    ylabel('SSRT - mean RT (ms)')
    
    
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
    
    %get inhibition function data with four values all together (maybe should
    %include those with three values as well)
    fourdigs=((arrayfun(@(x) length(x.inhibfun), alldata))==4 & (arrayfun(@(x) length(x.ssd), alldata))==4); %take only recordings where inhibition function and ssd have four values
    %threedigs=((arrayfun(@(x) length(x.inhibfun), alldata))==3 & (arrayfun(@(x) length(x.ssd), alldata))==3);
    fourdata=alldata(fourdigs);
    if ~isempty(fourdata)
        % figure;
        % plot([fourdata.ssd],[fourdata.inhibfun]);
        %same thing but rescaled SSD from 0 to 100
        fourdatassd=[fourdata.ssd];
        [~,sortdat]=sort(fourdatassd(4,:));
        % sfourdatassd=fourdatassd(:,sortdat);
        % for resc=1:size(sfourdatassd,2)
        %     sfourdatassd(:,resc)=sfourdatassd(:,resc)-sfourdatassd(1,resc);
        %     sfourdatassd(:,resc)=round(sfourdatassd(:,resc)./(sfourdatassd(4,resc)/100));
        % end
        fourdatainhibfun=[fourdata.inhibfun];
        sfourdatainhibfun=fourdatainhibfun(:,sortdat);
        % figure
        % subplot(2,1,1)
        % plot(cellfun(@(x) mean(x),fourdatanssdelay),'linewidth',2,'LineStyle','-.','color',[0.2 0.5 0.8])
        % hold on
        % plot(fourdatassd(4,:),'linewidth',2,'LineStyle','-.','color',[0.8 0.5 0.2])
        % subplot(2,1,2)
        
        %history of inhibition functions
        allmat=nan(400,size(fourdatassd,2)); % min(min(fourdatainhibfun)).*ones
        for session=1:size(allmat,2)
            allmat(fourdatassd(1,session)-50:fourdatassd(1,session),session)=fourdatainhibfun(1,session);
            allmat(fourdatassd(1,session)+1:fourdatassd(2,session),session)=fourdatainhibfun(2,session);
            allmat(fourdatassd(2,session)+1:fourdatassd(3,session),session)=fourdatainhibfun(3,session);
            allmat(fourdatassd(3,session)+1:fourdatassd(4,session),session)=fourdatainhibfun(4,session);
            %     allmat(fourdatassd(3,session):end,session)=fourdatainhibfun(4,session);
        end
        figure
        mesh(allmat)
        
        %sort nss sac delays
        fourdatanssdelay=cellfun(@(x) x.all, {fourdata(~cellfun('isempty',{fourdata.NSSsacdelay})).NSSsacdelay}, 'UniformOutput', false);
        sfourdatanssdelay=fourdatanssdelay(:,sortdat);
        
        %now split according to longest SSD
        sfourdatassd=fourdatassd(:,sortdat);
        binnum=sum(hist(fourdatassd(4,:))>1);
        if ~binnum
            binnum=sum(hist(fourdatassd(4,:))>0);
        end
        binsize=[0 hist(fourdatassd(4,:),binnum)];
        fdsplitif=cell(binnum,1);
        fdsplitssd=cell(binnum,1);
        fdsplitnsssacd=cell(binnum,1);
        
        % subplot to draw lines
        subplot(2,2,1);
        cc=lines(binnum);
        if size(cc,1)>=8
            cc(8,:)=[0 0.75 0];
        end
        hold on;
        
        for bins=2:binnum+1
            fdsplitif{bins-1}=sfourdatainhibfun(:,sum(binsize(1:bins-1))+1:sum(binsize(1:bins)));
            fdsplitssd{bins-1}=sfourdatassd(:,sum(binsize(1:bins-1))+1:sum(binsize(1:bins)));
            fdsplitnsssacd{bins-1}=[sfourdatanssdelay{:,sum(binsize(1:bins-1))+1:sum(binsize(1:bins))}];
            %         subplot(floor(binnum/2),2,bins-1)
            %         plot(fdsplitssd{bins-1},fdsplitif{bins-1});
        end
        SSDcats=vertcat(cellfun(@(x) [min(x(1,:)) max(x(4,:))], fdsplitssd,'UniformOutput', false));
        
        % now fit, and plot, each category with sigmoid
        yfitval=cell(binnum,1);
        hyfitvalline=nan(binnum,1);
        [mIntSSRT, mSSRT, oMeanSSRT]=deal(nan(binnum,2));
        for bins=1:binnum
            % plot mean of actual data
            splitssdbins=round(mean(fdsplitssd{bins},2));
            splitif=mean(fdsplitif{bins},2);
            plot(splitssdbins,splitif,'Color',cc(bins,:),'LineWidth',1.5);
            % calculate sigmoid fit
            fittime=0:10:500;
            fitresult = sigmoidfit(fdsplitssd{bins}, fdsplitif{bins});
            yfitval{bins}=fitresult(fittime); % This gives us the templates for the SSD range-dependant inhibition functions (six for Sixx, ha ha)
            % do some alignement between sigmoid and actual data
            ifstart=fittime(find(round(yfitval{bins}.*100)./100>=mean(fdsplitif{bins}(1,:),2),1)); %finds when the sig has same value has mean first inhibfun value
            ssdstart=fittime(find(fittime>mean(fdsplitssd{bins}(1,:),2),1));% fnd when mean first ssd value occurs
            sigdiff=ifstart-ssdstart;
            try
                fittime=fittime-sigdiff; %shift time
            catch
                sigdiff
            end
            % plot sigmoid
            hyfitvalline(bins)=plot(fittime(fittime>0),yfitval{bins}(fittime>0),'Color',cc(bins,:),'LineStyle','.','LineWidth',0.8);
            
            % and calculate SSD range-specific SSRTs
            [mIntSSRT(bins,1), mSSRT(bins,1), oMeanSSRT(bins,1)]= ...
                ssrt_bestfit(fdsplitnsssacd{bins}, splitif', splitssdbins');
            %             text(350,(binnum+1)/10-(bins/10),[num2str(round(mIntSSRT(bins,1))) ' '...
            %                 num2str(round(mSSRT(bins,1))) ' ' num2str(round(oMeanSSRT(bins,1)))],...
            %                 'Color',cc(bins,:),'FontName','calibri','FontSize',10);
            
            [mIntSSRT(bins,2), mSSRT(bins,2), oMeanSSRT(bins,2)]= ...
                ssrt_bestfit(fdsplitnsssacd{bins}, (yfitval{bins}(fittime>0))', fittime(fittime>0));
            %             text(480,(binnum+1)/10-(bins/10),[num2str(round(mIntSSRT(bins,2))) ' '...
            %                 num2str(round(mSSRT(bins,2))) ' ' num2str(round(oMeanSSRT(bins,2)))],...
            %                 'Color',cc(bins,:),'FontName','calibri','FontSize',8);
        end
        
        title('Sigmoid fit for each SSD range','FontName','calibri','FontSize',14);
        hxlabel=xlabel(gca,'SSDs','FontName','calibri','FontSize',12);
        set(gca,'Xlim',[50 400],'XTick',[100:50:350],'TickDir','out','box','off');
        hylabel=ylabel(gca,'P(Signal-respond)','FontName','calibri','FontSize',12);
        set(gca,'Ylim',[0 1],'TickDir','out','box','off');
        hlegif = legend(hyfitvalline, num2str(vertcat(SSDcats{:})),'Location','SouthEast');
        hlegifpos=get(hlegif,'position');
        hlegifpos(1)=0.5;
        set(hlegif,'Interpreter','none', 'Box', 'off','LineWidth',1.5,'FontName','calibri','FontSize',9,'position',hlegifpos);
        
        text(350,0.5,'data mean SSRT','Color','k','FontName','calibri','FontSize',10);
        text(350,0.4,[num2str(round(mean(mIntSSRT(:,1)))) ' '...
            num2str(round(mean(mSSRT(:,1)))) ' ' num2str(round(mean(oMeanSSRT(:,1))))],...
            'Color','k','FontName','calibri','FontSize',10);
        text(350,0.3,'sig fit mean SSRT', 'Color','r','FontName','calibri','FontSize',10);
        text(350,0.2,[num2str(round(mean(mIntSSRT(:,2)))) ' '...
            num2str(round(mean(mSSRT(:,2)))) ' ' num2str(round(mean(oMeanSSRT(:,2))))],...
            'Color','r','FontName','calibri','FontSize',10);
        
    end
    %plot meanssrt values
    %     figure;
    %     plot([1:length([allmeanssrt{~cellfun(@isempty,allmeanssrt)}])],[allmeanssrt{~cellfun(@isempty,allmeanssrt)}],'mo',...
    %     'LineWidth',1,...
    %     'MarkerEdgeColor','k',...
    %     'MarkerFaceColor',[1 .63 .49],...
    %     'MarkerSize',3);
    
    %% calculating probabilities
    
    allnccssd=vertcat(allnccssd{:});
    allccssd=vertcat(allccssd{:});
    allsacdelay=cellfun(@(x) x.all, {alldata(~cellfun('isempty',{alldata.NSSsacdelay})).NSSsacdelay}, 'UniformOutput', false); % concatenation super casse-c#*&!lle
    allsacdelay=[allsacdelay{:}];
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
    % fit sigmoid through it
    fitresult = sigmoidfit(narssdbins, probaresp);
    yfitval=fitresult(50:10:400); % This gives us the templates for the SSD range-dependant inhibition functions (six for Sixx, ha ha)
    
    % plot overall inhibition function
    subplot(2,2,3);
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
        ssrt_bestfit(allsacdelay, probaresp', narssdbins');
    text(0.5,0.6,['meanIntSSRT = ' num2str(round(meanIntSSRT))],'FontName','calibri','FontSize',10,'Parent',hvaldisp);
    text(0.5,0.5,['meanSSRT = ' num2str(round(meanSSRT))],'FontName','calibri','FontSize',10,'Parent',hvaldisp);
    text(0.5,0.4,['overallMeanSSRT = ' num2str(round(overallMeanSSRT))],'FontName','calibri','FontSize',10,'Parent',hvaldisp);
    % now with inhib function sigmoid fit
    [meanIntSSRT, meanSSRT, overallMeanSSRT]= ...
        ssrt_bestfit(allsacdelay, yfitval', [50:10:400]);
    text(0.5,0.3,['meanIntSSRT_s = ' num2str(round(meanIntSSRT))],'FontName','calibri','FontSize',10,'Parent',hvaldisp);
    text(0.5,0.2,['meanSSRT_s = ' num2str(round(meanSSRT))],'FontName','calibri','FontSize',10,'Parent',hvaldisp);
    text(0.5,0.1,['overallMeanSSRT_s = ' num2str(round(overallMeanSSRT))],'FontName','calibri','FontSize',10,'Parent',hvaldisp);
    
    
    %prepare and plot saccade latency frequency
    [saclatquant,saclatxlims]=hist(allsacdelay,[0:25:500]);
    saclatfreq=saclatquant./sum(saclatquant);
    % intsaclatfreq=interp(saclatfreq,2);
    % intsaclatxlims=interp(saclatxlims,2);
    subplot(2,2,4);
    plot(saclatxlims,saclatfreq,'Color','k','LineWidth',1.8);
    hold on
    %plot(intsaclatxlims,intsaclatfreq,'Color',[0.55 0.25 0.25],'LineWidth',1.8);
    if exist('fdsplitnsssacd') && size(fdsplitnsssacd,1)>1
        for bins=1:binnum
            [saclatquant,saclatxlims]=hist(fdsplitnsssacd{bins},[0:25:500]);
            saclatfreq=saclatquant./sum(saclatquant);
            %         intsaclatfreq=interp(saclatfreq,2);
            %         intsaclatxlims=interp(saclatxlims,2);
            plot(saclatxlims,saclatfreq,'Color',cc(bins,:),'LineWidth',1);
        end
    end
    title('Saccade Latency','FontName','calibri','FontSize',15);
    hxlabel=xlabel(gca,'Saccade latency','FontName','calibri','FontSize',12);
    set(gca,'Xlim',[0 500],'XTick',[0:50:500],'TickDir','out','box','off'); %'XTickLabel',[50:50:400]
    hylabel=ylabel(gca,'Proportion','FontName','calibri','FontSize',12);
    curylim=get(gca,'YLim');
    set(gca,'Ylim',[0 curylim(2)],'TickDir','out','box','off');
    text(350,curylim(2)/1.5,'NSS success rate','FontName','calibri','FontSize',10);
    text(400,curylim(2)/1.8,num2str(mean([alldata.NSSsuccessrate])),'FontName','calibri','FontSize',10);
    exportfigname=['CMdat_',subject,'_',emdirections{andir}];
    print(gcf, '-dpng', '-noui', '-opengl','-r600', exportfigname);
    delete(CMdatfig);
    
    % figure
    % plot(narssdbins,probaresp,'o','MarkerEdgeColor','k','MarkerFaceColor','c','MarkerSize',8)
    % set(gca,'XTick',[1:24],'XTickLabel',narssdbins)
    % Generalized linear model regression
    % regcoeffs = glmfit(narssdbins,[nccssdhist' (nccssdhist'+ccssdhist')],'binomial','link','probit');
    % respregfit = glmval(regcoeffs, narssdbins,'probit','size', (nccssdhist'+ccssdhist'));
    % figure
    % subplot(2,2,3)
    % plot(narssdbins, nccssdhist'./(nccssdhist'+ccssdhist'),'o',...
    %     narssdbins,respregfit./(nccssdhist'+ccssdhist'),'-','LineWidth',2)
    % xlabel('stop signal delay (ms)');
    % title('Probability to respond')
    %
    % %% finding ssrt - first method: integration (constant SSRT)
    % delaydistribedges=min(allsacdelay)-1:10:max(allsacdelay)+1;
    % delaydistribhhist=histc(sort(allsacdelay),delaydistribedges);
    % emptyvalues=find(~delaydistribhhist);
    % if emptyvalues(1)==1
    %     emptyvalues=emptyvalues(2:end);
    % end
    % if emptyvalues(end)==length(delaydistribhhist)
    %     emptyvalues=emptyvalues(1:end-1);
    % end
    % for empt=1:length(emptyvalues)
    %     delaydistribhhist(emptyvalues(empt))=mean([delaydistribhhist(emptyvalues(empt)-1)...
    %         delaydistribhhist(emptyvalues(empt)+1)]);
    % end
    % %interpdelaydistribpdf=interp1(delaydistribedges,delaydistribpdf,min(totsadelay)-1:1:max(totsadelay)+1);
    % delayfreq=delaydistribhhist./sum(delaydistribhhist);
    % subplot(2,2,1)
    % plot(delaydistribedges,delayfreq,'-b','LineWidth',2);
    % xlabel('Saccade delay (ms)');
    % title('Saccade delay frequency for no-stop trials')
    % subplot(2,2,2)
    % plot(delaydistribedges,cumtrapz(delayfreq),'-','LineWidth',2);
    % xlabel('Saccade delay (ms)');
    % title('Cumulative distribution of saccade delay')
    %
    % %interpdelfrq=interpdelaydistribpdf./sum(interpdelaydistribpdf);
    % %     wbfitparam=wblfit(delayfreq);
    % %     delaywbfit=wblcdf(delayfreq,wbfitparam(1),wbfitparam(2));
    %
    % % finding SSRT for each SSD
    % intssrt=nan(length(narssdbins),1);
    % for prbr=1:length(narssdbins)
    % sufficientdelay=find(probaresp(prbr)<=cumtrapz(delayfreq),1);
    % stopprocfinishline=delaydistribedges(sufficientdelay);
    % intssrt(prbr)=stopprocfinishline-narssdbins(prbr);
    % end
    %
    %
    % %% finding ssrt - second method (random SSRT)
    %
    % %mean inhibition function
    % mininhibfun=(sum((probaresp(2:end)-probaresp(1:end-1)).*narssdbins(2:end)))./(max(probaresp)-min(probaresp));
    % inhibfunssrt(1)=mean(allsacdelay)-mininhibfun;
    % % fitting weibull function
    % if ~(length(find(probaresp))<length(probaresp))
    %     weibullcdfparam=wblfit(probaresp);
    %     weibullcdf=wblcdf(probaresp,weibullcdfparam(1),weibullcdfparam(2));
    %     %plot(weibullcdf)
    %     mininhibfunwb=(sum((weibullcdf(2:end)-weibullcdf(1:end-1)).*narssdbins(2:end)))...
    %         ./(max(weibullcdf)-min(weibullcdf));
    %     inhibfunssrt(2)=mean(allsacdelay)-mininhibfunwb;
    % end
    %
    % %% averaging different estimates
    % if mean(inhibfunssrt)>=(mean(intssrt)+std(intssrt)) || mean(inhibfunssrt)<=(mean(intssrt)-std(intssrt))
    %     ssrt=mean(intssrt);
    % else
    %     %ssrt=mean([intssrt' inhibfunssrt]);
    %     ssrt=mean([mean(intssrt) inhibfunssrt]);
    % end
    %
    % subplot(2,2,4)
    % plot(narssdbins,intssrt,'-','LineWidth',2)
    % hold on
    % plot(narssdbins,ones(length(narssdbins),1)*mean(inhibfunssrt),'-g','LineWidth',2)
    % plot(narssdbins,ones(length(narssdbins),1)*ssrt,'-r','LineWidth',2)
    % legend({'constant SSRT method','random SSRT method','average different estimates'})
    % xlabel('stop signal delay (ms)');
    % title('Stop signal reaction time');
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
    
    %% output file
    
    % txtfilename = ['SSRTs_',subject,'.txt'];
    % fid = fopen(txtfilename, 'a');
    % fprintf(fid,'%10s\t %10s\t %10s\t %10s\t %10s\r\n','latencies', 'durations', 'amplitudes', 'peakvels', 'gapdelays');
    % fprintf(fid, '%10.0f\t %10.0f\t %10.1f\t %10.0f\t %10.0f\r\n', [latencies'; durations'; amplitude'; pkvels'; gapdelays']);
    % fclose(fid);
    
    %% collect alldata
    if ~exist('colecalldata')
        %colecalldata=struct('dir',[],'NSSsacdelay',[],'NSSsuccessrate',[],'SRsacdelay',[],'ssd',[],'inhibfun',[],'meanIntSSRT',[],'meanSSRT',[],'overallMeanSSRT',[]);
        colecalldata=alldata;
    else
        colecalldata=[colecalldata alldata];
    end
end
end


