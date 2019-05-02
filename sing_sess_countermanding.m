function sing_sess_countermanding(fileInfo,cmdData)
% Analyze and display data from a single session of countermanding task

global directory slash;

%% settings
userinfo=SetUserDir;
directory=userinfo.directory;
slash=userinfo.slash;
% CmdData=struct('file',[],'ssrt',[],'tach',[],'saclat',[],'sdf',[],...
%     'rast',[],'algidx',[],'allviscuetimes',[]);
if ~exist('cmdData','var')
    cmdData=struct('file',[],'allalignmnt',[],'allprevssd',[],'allssds',[],'allsacdelay',[],...
        'allprefdir',[],'allndata',[],'allmssrt_tacho',[],'alldb',[],'normFactor',[]);
    cd(userinfo.syncdir);
end
if strcmp(fileInfo{1},'get')
    %% finding best candidate
    
    % userinfo.user,userinfo.dbldir,userinfo.mapdr,userinfo.servrep,userinfo.mapddataf
    CCNdb = connect2DB('vp_sldata');
    
    cd(userinfo.syncdir);
    if ~exist('cmdData','var')
        load('cDn_gsdata.mat'); %cDn_gsdata.mat  top_cortex_gsdata.mat
    end
    %number cells
    goodrecs=~cellfun('isempty',cmdData.allsacdelay);
    % st.goodrecs=~cellfun('isempty',stdata.allsacdelay);
    
    %remove bad apples
    fieldName = fieldnames(cmdData);
    for lp=1:length(fieldName)
        cmdData.(fieldName{lp})=cmdData.(fieldName{lp})(goodrecs,:);
    end
    
    sigma=15;
    periEventActivity=cellfun(@(x) conv_raster(x(1,1).rast,sigma,x(1,1).alignt-(300+sigma*3),x(1,1).alignt+(299+sigma*3)), cmdData.allndata(:,1), 'UniformOutput',false); %600ms epoch
    
    % remove recordings with one trial
    properLengthRecs=cellfun(@(x) length(x)>1,periEventActivity);
    goodrecs(goodrecs)=properLengthRecs;
    periEventActivity=periEventActivity(properLengthRecs);
    
    % get unit cluster info and profiles
    unit_ids=cellfun(@(x) x.unit_id,cmdData.alldb);
    
    [sorted_unit_ids,sunitid_idx]=sort(unit_ids);
    query = ['SELECT c.profile, c.profile_type FROM clusters c WHERE cluster_id IN (' sprintf('%.0f,' ,sorted_unit_ids(1:end-1)') num2str(sorted_unit_ids(end)) ')'];
    profiles = fetch(CCNdb,query);
    sunitid_revidx(sunitid_idx)=1:length(unit_ids);
    clusidx=[profiles{sunitid_revidx,2}];
    clustypes={profiles{sunitid_revidx,1}};
    
    clus1Idx=clusidx==101;
    
    % find most active cells, aligned to saccade, in cluster #1
    periEventActivity=periEventActivity(clus1Idx);
    maxEvent=nan(size(periEventActivity,1),1);
    maxEvent=cellfun(@(eventEpoch) max(eventEpoch),periEventActivity, 'UniformOutput',true);
    
    [~,topEventActivity]=sort(maxEvent,'descend');
    
    % figure;hold on
    % for ev=1:10
    %     plot(periEventActivity{topEventActivity(ev)});
    % end
    % legend(num2str((1:10)'))
    
    [~,topEventActivity]=sort(topEventActivity);
    
    % hand picked most active cells: topEventActivity([2,3,5,8,10]) ->     14,5,4,31,15
    
    % now gets those with the most trials, aligned to target
    numTrials=cellfun(@(x) [size(x(1).rast,1) size(x(2).rast,1) size(x(3).rast,1)],cmdData.allndata(:,2), 'UniformOutput',false);
    % mostTrials=cat(1,numTrials{[14,5,4,31,15]});
    numTrials=cat(1,numTrials{:});
    numTrials=numTrials(clus1Idx,:);
    [~,mostNSST]=sort(numTrials(:,1),'descend'); [~,mostNSSTscore]=sort(mostNSST);
    [~,mostCSST]=sort(numTrials(:,2),'descend'); [~,mostCSSTscore]=sort(mostCSST);
    [~,mostNCSST]=sort(numTrials(:,3),'descend'); [~,mostNCSSTscore]=sort(mostNCSST);
    %     Index can do without mostNSSTscore
    [~,mostTrials]=sort(sum([mostNSSTscore,mostCSSTscore,mostNCSSTscore],2));[~,mostTrials]=sort(mostTrials);
    %     [~,mostTrials]=sort(sum([mostCSSTscore,mostNCSSTscore],2));[~,mostTrials]=sort(mostTrials);
    
    % make an index of recordings with most same-ssd in SST.
    % now gets those with the most trials, aligned to target
    numSSDs=cellfun(@(x) [x{1,1};x{1,2}] ,cmdData.allssds, 'UniformOutput',false);
    numSSDs=numSSDs(clus1Idx,:);
    [~,numSSDs]=sort(cellfun(@(x) length(unique(floor(x/5)*5)), numSSDs, 'UniformOutput',true));[~,numSSDs]=sort(numSSDs);
    
    %combine indices to find recordings with best activity and most trials
    [~,bestFiles]=sort(topEventActivity+mostTrials+numSSDs);
    %     numTrials(bestFiles,:)
    % hand picked recordings with best activity and most trials were 5,14,15
    
    % plot them
    %     for cellNum=1:10
    %         bestFile=bestFiles(cellNum);
    %
    %         figure('position',[1686 49 867 1308])
    %         for plotNum=1:3
    %             %get rasters
    %             periTargetActivity=cellfun(@(x) conv_raster(x(plotNum).rast,sigma,x(plotNum).alignt-(100+sigma*3),x(plotNum).alignt+(500+sigma*3)), cmdData.allndata(:,2), 'UniformOutput',false); %600ms epoch
    %             periTargetActivity=periTargetActivity(clus1Idx);
    %
    %             subplot(3,1,plotNum)
    %             plot(periTargetActivity{bestFile},'LineWidth',2);
    %             currylim=get(gca,'YLim');
    %                 patch([repmat(98,1,2) repmat(102,1,2)], ...
    %                 [[0 currylim(2)] fliplr([0 currylim(2)])], ...
    %                 [0 0 0 0],[1 0 0],'EdgeColor','none','FaceAlpha',0.5);
    %
    %         axis(gca,'tight');
    %         box off;
    %
    %         end
    %     end
    
    %best cells in cDN cluster #1 are 31, 4, 8, 1
    
    clus1Idx=find(clus1Idx);
    fieldName = fieldnames(cmdData);
    for lp=1:length(fieldName)
        cmdData.(fieldName{lp})=cmdData.(fieldName{lp})(clus1Idx(bestFiles(1)),:);
    end
    % queries{1} = ['SELECT a_file FROM sorts s INNER JOIN recordings r on s.recording_id_fk = r.recording_id WHERE sort_id IN (' ...
    %     sprintf('%.0f,' ,cellfun(@(x) x.sort_id,cmdData.alldb(1:end-1,1))) num2str(cmdData.alldb{end,1}.sort_id) ')'];
    cmdData.fileName=fetch(CCNdb,['SELECT a_file FROM sorts s INNER JOIN recordings r on ' ...
        's.recording_id_fk = r.recording_id WHERE sort_id IN (' ...
        num2str(cmdData.alldb{clus1Idx(bestFiles(1))}.sort_id) ')']);
    cmdData.fileName=cmdData.fileName{:}(1:end-1);
    
else
    cmdData.fileName=fileInfo{1};
end

if strcmp('R',cmdData.fileName(1))
    subject='Rigel';
    procdir = [directory,'processed',slash,'Rigel',slash];
elseif strcmp('S',cmdData.fileName(1))
    subject='Sixx';
    procdir = [directory,'processed',slash,'Sixx',slash];
elseif strcmp('H',cmdData.fileName(1))
    subject='Hilda';
    procdir = [directory,'processed',slash,'Hilda',slash];
end

procdirlisting=dir(procdir);
procdirfileNames={procdirlisting.name};
loadfile=procdirfileNames(~cellfun('isempty',regexpi(procdirfileNames,cmdData.fileName,'match')));
% if two versions (REX and Spike2), chose Spike2
if length(loadfile)>1
    loadfile=loadfile(~cellfun('isempty',regexpi(loadfile,'Sp2','match')));
    if isempty(loadfile) %actually it's an old file that sneaked in
        loadfile=procdirfileNames(~cellfun('isempty',regexpi(procdirfileNames,cmdData.fileName,'match')));
        loadfile=loadfile(~cellfun('isempty',regexpi(loadfile,'REX','match')));
    end
end

if isempty(cmdData.allndata) %redundant now
    cmdData=load('cDn_gsdata.mat','gsdata'); %cDn_gsdata.mat  top_cortex_gsdata.mat
    cmdData=cmdData.gsdata;
end
if size(cmdData.allndata,2)>1
    cellIdx=ismember(cmdData.alldb.rec_id,fileInfo{2}); %_id cellfun(@(dbInfo) dbInfo.rec_id==fileInfo{2}, cmdData.alldb); %str2double
    fn = fieldnames(cmdData); fn=fn(~cellfun(@(fieldname) strcmp(fieldname,'fileName'),fn));
    for lp=1:length(fn)
        try
            cmdData.(fn{lp})=cmdData.(fn{lp})(cellIdx,:);
        catch
            [cmdData.(fn{lp})]=deal(cmdData.(fn{lp})(cellIdx));
        end
    end
    %     cmdData.fileName=fileInfo{1};
end

if isempty(cmdData.allmssrt_tacho)
    %% get countermanding session results
    [mssrt,inhibfun,ccssd,nccssd,ssdvalues,tachomc,tachowidth,sacdelay,rewtimes]=findssrt(loadfile{:}, 0);
end

%% session psychophysics
if ~isfield(cmdData,'allssds')
    
    % subjects={'Rigel','Sixx','Hilda'};
    %
    % % best use a recording made with gapstop training and lots of trials
    % subject=subjects{3};
    % recname='H48L5A1_18102'; % S148cnttrain
    % load(recname,'allcodes','alltimes','allbad','saccadeInfo');
    % % load([recname,'_sac'],'alignedata');
    %
    % %% find stop trials
    %     trialtypes=floor(allcodes(:,2)./10);%./10 only /10 if pooling data together
    %     stoptrials=find(trialtypes==407);
    %     stoptrialcodes=allcodes(stoptrials,:);
    %
    % %% find canceled and noncanceled stop trials
    % if find(stoptrialcodes(:,8)==1503,1) %recordings with benchmark code (1503)
    %     bmssd=1;
    %     abortedstoptrials=~(floor(stoptrialcodes(:,9)./10)==507);
    %     noncancel=(stoptrialcodes(~abortedstoptrials,10)==16386) | ... % trials where a saccade is initiated
    %         (stoptrialcodes(~abortedstoptrials,10)==17385);                % or subsequently breaking fixation
    % else
    %     bmssd=0;
    %     abortedstoptrials=~(floor(stoptrialcodes(:,8)./10)==507);
    %     noncancel=(stoptrialcodes(~abortedstoptrials,9)==17385) | ...
    %         (stoptrialcodes(~abortedstoptrials,9)==16386);
    % end
    %
    % %% find "desired" delay times (keep in mind that video synch adds ~65ms, so real ssd are variable)
    %
    % % calculate ssd
    % stoptrialtimes=alltimes(stoptrials(~abortedstoptrials,:),:);
    % noncanceltimes=stoptrialtimes(noncancel,:);
    % canceltimes=stoptrialtimes(~noncancel,:);
    %
    %     if bmssd %benchmark
    %         nccssd=(noncanceltimes(:,9)-noncanceltimes(:,7))-3; %3 ms added by the state transitions
    %         ccssd=(canceltimes(:,9)-canceltimes(:,7))-3;
    %         allssd=nan(size(alltimes,1),1);
    %         allssd(stoptrials)=(alltimes(stoptrials,9)-alltimes(stoptrials,7))-3;
    %     else
    %         nccssd=(noncanceltimes(:,8)-noncanceltimes(:,7))-2; %2 ms added by the state transitions
    %         ccssd=(canceltimes(:,8)-canceltimes(:,7))-2;
    %     end
    %
    % ssdvalues=unique([ccssd;nccssd]);
    %
    % % show progression of rt and ssd
    % % figure;
    % % subplot(2,1,1);
    % % stairs(goodtrials);
    % % set(gca,'xlim',[1 length(allssd)]);
    % % title('successful trials');
    % % allssd=zeros(size(alltimes,1),1);
    % % allssd(stoptrials(~abortedstoptrials))=...
    % %     (alltimes(stoptrials(~abortedstoptrials),9)-alltimes(stoptrials(~abortedstoptrials),7))-3;
    % % subplot(2,1,2);
    % % plot(find(allssd>0),allssd(allssd>0));
    % % set(gca,'xlim',[1 length(allssd)]);
    % % hold on;
    %
    % %% get saccade delay for signal-respond (noncancel) trial
    %
    %     if bmssd %benchmark
    %         alldata.SRsacdelay=(noncanceltimes(:,10)-noncanceltimes(:,7))-6; %6 ms added by the state transitions
    %     else
    %         alldata.SRsacdelay=(noncanceltimes(:,9)-noncanceltimes(:,7))-6;
    %     end
    %
    % %% saccade delay for non-stop trials: all good saccade from non-stop trials
    % %(may yield slightly different results than with left/right parsing method
    % % used previously)
    %
    % alllats=reshape({saccadeInfo.latency},size(saccadeInfo));
    % alllats=alllats';%needs to be transposed because the logical indexing below will be done column by column, not row by row
    % allgoodsacs=~cellfun('isempty',reshape({saccadeInfo.latency},size(saccadeInfo)));
    %     %removing bad trials
    %     allgoodsacs(logical(allbad),:)=0;
    %     %removing stop trials that may be included
    %     allgoodsacs(floor(allcodes(:,2)./1000)~=6,:)=0;
    %     %indexing good sac trials
    %     % if saccade detection corrected, there may two 'good' saccades
    %     if max(sum(allgoodsacs,2))>1
    %         twogoods=find(sum(allgoodsacs,2)>1);
    %         for dblsac=1:length(twogoods)
    %             allgoodsacs(twogoods(dblsac),find(allgoodsacs(twogoods(dblsac),:),1))=0;
    %         end
    %     end
    % sacdelay.all=(cell2mat(alllats(allgoodsacs')))';
    %
    % alldata.NSSsacdelay=sacdelay;%[sacdelay{:}];
    %
    % % plot RTs on top of SSDs
    % % alllats=zeros(size(allgoodsacs,1),1);
    % % alllats(logical(sum(allgoodsacs,2)))=sacdelay.all;
    % % plot(find(alllats>0),alllats(alllats>0),'r');
    % % title('evolution of SSDs and RTs across trials')
    % % legend('ssd','rt');
    %
    % % NSS success rate
    % % NSS_targ=allcodes(floor(allcodes(:,7)./10)==684,8);
    % % alldata.NSSsuccessrate=sum(floor(NSS_targ./10)==704)/length(NSS_targ);
    %
    % %% calculating probabilities for this session
    % %histo binning method
    % [~,delaybincenters]=hist([ccssd;nccssd],4);
    % alldlbincnt=round(delaybincenters');
    % % sort unique delay bin centers
    % ssdbins=unique(sort(alldlbincnt));
    % % narrow ssds to those found more than once
    %     % narssdbins=ssdbins(hist(nccssd,ssdbins)>1); %not for
    %     % individual SSRT calculation
    %     narssdbins=ssdbins(hist(nccssd,ssdbins)>0);
    %     alldata.ssd=narssdbins;
    % % bin canceled and non-canceled trials according to these ssds
    % try
    % nccssdhist=hist(nccssd,narssdbins);
    % ccssdhist=hist(ccssd,narssdbins);
    % % find probability to respond
    % probaresp=nccssdhist'./(nccssdhist'+ccssdhist');
    % catch
    %     probaresp=1;
    %     narssdbins=0;
    % end
    %
    % %diff method
    % ssdvalues(find(diff(ssdvalues)==1)+1)=ssdvalues(diff(ssdvalues)==1);
    % ssdvalues=ssdvalues(diff(ssdvalues)>0);
    % try
    % nccssdhist=hist(nccssd,ssdvalues);
    % ccssdhist=hist(ccssd,ssdvalues);
    % % find probability to respond
    % probaresp_diff=nccssdhist'./(nccssdhist'+ccssdhist');
    % catch
    %     probaresp_diff=1;
    % end
    %
    % %% calculate SSRT
    % if ~(isempty(narssdbins) || length(narssdbins)==1)
    %
    %     % test monotonicity and keep relevant inhibition function
    % if (all(diff(probaresp)>=0) && length(probaresp)>=3) && ~all(diff(probaresp_diff)>=0)
    %     alldata.inhibfun=probaresp;
    %         % calculate SSRT with Boucher et al's method
    %         try
    %             [meanIntSSRT, meanSSRT, overallMeanSSRT]= ...
    %             ssrt_bestfit(sacdelay.all', probaresp', narssdbins');
    %         catch
    %            sacdelay;
    %         end
    % else
    %     alldata.inhibfun=probaresp_diff;
    %         % calculate SSRTs
    %         try
    %         [meanIntSSRT, meanSSRT, overallMeanSSRT]= ...
    %             ssrt_bestfit(sacdelay.all', probaresp_diff', ssdvalues');
    %         catch
    %         end
    % end
    %
    % allmeanssrt=overallMeanSSRT;
    % alldata.meanIntSSRT=meanIntSSRT;
    % alldata.meanSSRT=meanSSRT;
    % alldata.overallMeanSSRT=overallMeanSSRT;
    % end
    % % clearvars -except nccssd ccssd allsacdelay alldlbincnt filestoload ...
    % %     numfile splitdataaligned allmeanssrt alldata andir emdirections subject ...
    % %     monknum colecalldata;
    %
    % %% prepare plot
    %     CMdatfig=figure;
    % %     CMdatfigpos=get(CMdatfig,'Position');
    %     CMdatfigpos=[500 300 800 500];
    %     set(CMdatfig,'Position',CMdatfigpos,'Color','w');
    %
    %         % subplot to display values
    %     hvaldisp = axes('Position', [.6, .6, .3, .3]);
    %     set(hvaldisp,'Visible','off');
    %         htitle=title('SSD ranges and SSRT values','FontName','calibri','FontSize',14);
    %         set(htitle,'Visible','on')
    %
    % % fit sigmoid through inhibition function
    %         fitresult = sigmoidfit(narssdbins, probaresp);
    %         yfitval=fitresult(50:10:400); % This gives us the templates for the SSD range-dependant inhibition functions (six for Sixx, ha ha)
    %
    % % plot overall inhibition function
    % subplot(1,2,1);
    % plot(narssdbins,probaresp,'Color',[0.25 0.25 0.25],'LineWidth',1.8);
    % % hold on
    % % plot([50:10:400],yfitval,'Color',[0.2 0.4 0.6],'LineStyle','-.','LineWidth',1.5);
    % title('Inhibition function','FontName','calibri','FontSize',15);
    % hxlabel=xlabel(gca,'Stop Signal Delays','FontName','calibri','FontSize',12);
    % set(gca,'Xlim',[50 400],'XTick',[100:50:350],'TickDir','out','box','off'); %'XTickLabel',[50:50:400]
    % hylabel=ylabel(gca,'P(Signal-respond)','FontName','calibri','FontSize',12);
    % set(gca,'Ylim',[0 1],'TickDir','out','box','off');
    %
    % % print SSRT values
    % [meanIntSSRT, meanSSRT, overallMeanSSRT]= ...
    %     ssrt_bestfit(sacdelay.all', probaresp', narssdbins');
    % text(0.5,0.6,['meanIntSSRT = ' num2str(round(meanIntSSRT))],'FontName','calibri','FontSize',10,'Parent',hvaldisp);
    % text(0.5,0.5,['meanSSRT = ' num2str(round(meanSSRT))],'FontName','calibri','FontSize',10,'Parent',hvaldisp);
    % text(0.5,0.4,['overallMeanSSRT = ' num2str(round(overallMeanSSRT))],'FontName','calibri','FontSize',10,'Parent',hvaldisp);
    %     % now with inhib function sigmoid fit
    %     [meanIntSSRT, meanSSRT, overallMeanSSRT]= ...
    %     ssrt_bestfit(sacdelay.all', yfitval', [50:10:400]);
    %     text(0.5,0.3,['meanIntSSRT_s = ' num2str(round(meanIntSSRT))],'FontName','calibri','FontSize',10,'Parent',hvaldisp);
    %     text(0.5,0.2,['meanSSRT_s = ' num2str(round(meanSSRT))],'FontName','calibri','FontSize',10,'Parent',hvaldisp);
    %     text(0.5,0.1,['overallMeanSSRT_s = ' num2str(round(overallMeanSSRT))],'FontName','calibri','FontSize',10,'Parent',hvaldisp);
    %
    %
    % %prepare and plot saccade latency frequency
    % [saclatquant,saclatxlims]=hist(sacdelay.all,[0:25:500]);
    % saclatfreq=saclatquant./sum(saclatquant);
    %
    % subplot(1,2,2);
    % plot(saclatxlims,saclatfreq,'Color','k','LineWidth',1.8);
    % title('Saccade Latency','FontName','calibri','FontSize',15);
    % hxlabel=xlabel(gca,'Saccade latency','FontName','calibri','FontSize',12);
    % set(gca,'Xlim',[0 500],'XTick',[0:50:500],'TickDir','out','box','off'); %'XTickLabel',[50:50:400]
    % hylabel=ylabel(gca,'Proportion','FontName','calibri','FontSize',12);
    % curylim=get(gca,'YLim');
    % set(gca,'Ylim',[0 curylim(2)],'TickDir','out','box','off');
    %
    % %% export figure
    % exportfigname=['CMdat_',subject,'_',recname];
    % print(gcf, '-dpng', '-noui', '-opengl','-r600', exportfigname);
    % delete(CMdatfig);
end

%% session recording
singleSSD=true;
alignCondition='basic_tgt';

alignedata=cmdData.allndata.target; %{1, 2};
if strcmp(alignCondition,'basic_tgt')
    alignedata=alignedata(1:2);
end

numrast=length(alignedata);
[alignedata(1:numrast).alignlabel]=deal(cmdData.allalignmnt{1, 2});
% rasters=cmdData.allndata{1, 2}.rast;  %alignedata(rastnum).rasters;
% alignidx=cmdData.allndata{1, 2}.alignt;  %alignedata(rastnum).alignt;
% greyareas=cmdData.allndata{1, 2}.evttime; %alignedata(rastnum).allgreyareas

% for single ssd plots
if singleSSD
    
    % find most prevalent SSD
    %     [ssdbin,ssdbinval]=hist([cmdData.allssds{1, 1}{1, 1};cmdData.allssds{1, 1}{1, 2}]);
    %     ssdspread=abs(cmdData.allprevssd{1, 1}{:}-max(ssdbinval(ssdbin==max(ssdbin))));
    %     prevssd=cmdData.allprevssd{1, 1}{:}(ssdspread==min(ssdspread));
    allSSDs=[cmdData.allssds{1, 1}{1, 1};cmdData.allssds{1, 1}{1, 2}];
    [~,~,uniqSSDsIdx]=unique(floor(allSSDs/5)*5);
    % if plotting second most frequent
    %     prevssdIdx=uniqSSDsIdx==mode(uniqSSDsIdx);
    %     allSSDs=allSSDs(~prevssdIdx);uniqSSDsIdx=uniqSSDsIdx(~prevssdIdx);
    prevssd=mode(allSSDs(uniqSSDsIdx==mode(uniqSSDsIdx)));
    %categorize prevalent SSD
    SSDplace=prevssd/median(allSSDs);
    if SSDplace<0.7
        SSDcat='low range';
    elseif SSDplace>=0.7 & SSDplace<1.3
        SSDcat='mid range';
    elseif SSDplace>=1.3
        SSDcat='high range';
    end
    condition=['Single SSD (' num2str(prevssd) ') ' SSDcat];
    
    % CSST Latency matched NSS trials
    % Keeping NSS trials with sac latencies long enough
    % that they would have occured after the stop-signal and ssrt
    
    fieldName = fieldnames(alignedata);
    fieldName=fieldName(structfun(@(x) numel(x)>1 & ~ischar(x), alignedata(1, 1),'UniformOutput', true));
    
    for lp=1:length(fieldName)
        try
            alignedata(1, 1).(fieldName{lp})=alignedata(1, 1).(fieldName{lp})(...
                cmdData.allsacdelay{1, 1}.nsst>prevssd+round(cmdData.allmssrt_tacho{1, 1}{1}),:);
        catch
            alignedata(1, 1).(fieldName{lp})=alignedata(1, 1).(fieldName{lp})(...
                :,cmdData.allsacdelay{1, 1}.nsst>prevssd+round(cmdData.allmssrt_tacho{1, 1}{1}));
        end
    end
    
    %keep relevant SS trials (with SSD within +-3ms of prevalent SSD)
    %CS trials
    for lp=1:length(fieldName)
        try
            alignedata(1, 2).(fieldName{lp})=alignedata(1, 2).(fieldName{lp})(logical(arrayfun(@(x) sum(prevssd<=x+3 & prevssd>=x-3),...
                cmdData.allssds{1, 1}{1, 1})),:);
        catch
            alignedata(1, 2).(fieldName{lp})=alignedata(1, 2).(fieldName{lp})(:,logical(arrayfun(@(x) sum(prevssd<=x+3 & prevssd>=x-3),...
                cmdData.allssds{1, 1}{1, 1})));
        end
    end
    
    if ~strcmp(alignCondition,'basic_tgt')
        %NCS trials
        for lp=1:length(fieldName)
            try
                alignedata(1, 3).(fieldName{lp})=alignedata(1, 3).(fieldName{lp})(logical(arrayfun(@(x) sum(prevssd<=x+3 & prevssd>=x-3),...
                    cmdData.allssds{1, 1}{1, 2})),:);
            catch
                alignedata(1, 3).(fieldName{lp})=alignedata(1, 3).(fieldName{lp})(:,logical(arrayfun(@(x) sum(prevssd<=x+3 & prevssd>=x-3),...
                    cmdData.allssds{1, 1}{1, 2})));
            end
        end
    end
    %     alignedata(1, 3).rast=alignedata(1, 3).rast(logical(arrayfun(@(x) sum(prevssd<=x+3 & prevssd>=x-3),...
    %         cmdData.allssds{1, 1}{1, 2})),:);
else
    condition='All SSDs';
end

%plot related variables
plotstart=200;
plotstop=600;
fsigma=30;
causker=0;
chrono=1;

mainfig=figure('position',[1186 321 374 601]); cmap=colormap('lines');

for rastnum=1:numrast
    rasters=alignedata(rastnum).rast;
    alignidx=alignedata(rastnum).alignt;
    if chrono
        cut_chrasters=zeros(size(rasters,1),plotstart+plotstop+1);
        %listing relevant trials in a continuous series with other rasters
        chronoidx=ismember(sort([alignedata.trialnb]),alignedata(rastnum).trialnb);
    end
    greyareas=alignedata(rastnum).evttime;
    
    start=alignidx - plotstart;
    stop=alignidx + plotstop;
    
    if start < 1
        start = 1;
    end
    if stop > length(rasters)
        stop = length(rasters);
    end
    
    if chrono
        onerastplot=subplot(2,1,1,'Layer','top', ...
            'XTick',[],'YTick',[],'XColor','white','YColor','white', 'Parent', mainfig);
    else
        if numrast==1
            hrastplot(rastnum)=subplot(2,1,1:2,'Layer','top', ...
                'XTick',[],'YTick',[],'XColor','white','YColor','white', 'Parent', mainfig);
        else
            hrastplot(rastnum)=subplot(2,1,rastnum,'Layer','top', ...
                'XTick',[],'YTick',[],'XColor','white','YColor','white', 'Parent', mainfig);
        end
    end
    %reducing spacing between rasters
    if numrast>1 && ~chrono
        rastpos=get(gca,'position');
        rastpos(2)=rastpos(2)+rastpos(4)*0.5;
        set(gca,'position',rastpos);
    end
    
    if ~chrono
        % sorting rasters according greytime
        viscuetimes=nan(size(greyareas,2),2);
        for grst=1:size(greyareas,2)
            viscuetimes(grst,:)=greyareas{grst}(1,:);
        end
        cuestarts=viscuetimes(:,1);
        [~,sortidx]=sort(cuestarts,'descend');
        viscuetimes=viscuetimes(sortidx,:);
        rasters=rasters(sortidx,:);
    end
    
    hold on
    
    cut_rasters = rasters(:,start:stop); % Isolate rasters of interest
    cut_rast_siz = size(cut_rasters);
    isnantrial = isnan(sum(cut_rasters,2)); % Identify nantrials
    cut_rasters(isnan(cut_rasters)) = 0; % take nans out so they don't get plotted
    if chrono
        cut_chrasters(chronoidx,:)=cut_rasters;
        [indy, indx] = ind2sub(size(cut_chrasters),find(cut_chrasters)); %find row and column coordinates of spikes
    else
        [indy, indx] = ind2sub(size(cut_rasters),find(cut_rasters)); %find row and column coordinates of spikes
    end
    
    
    if(size(rasters,1) == 1) %don't plot
        %         plot([indx;indx],[indy;indy+1],'color',cmap(rastnum,:),'LineStyle','-'); % plot rasters
    else
        plot([indx';indx'],[indy';indy'+1],'color',cmap(rastnum,:),'LineStyle','-'); % plot rasters
    end
    
    % drawing the grey areas
    try
        greytimes=viscuetimes-start;
        greytimes(greytimes<0)=0;
        greytimes(greytimes>(stop-start+1))=stop-start+1;
    catch
        greytimes=0;
    end
    
    if ~sum(sum(isnan(greytimes))) && logical(sum(sum(greytimes))) && ~chrono
        grxlims=[greytimes';greytimes(:,2:-1:1)'];
        grylims=[1:size(grxlims,2);1:size(grxlims,2);2:size(grxlims,2)+1;2:size(grxlims,2)+1];
        patch(grxlims, grylims, [0 0 0], 'EdgeColor', 'none','FaceAlpha',0.2); % target grey area for targte aligment is not usefull
    end
    
    set(gca,'xlim',[1 length(start:stop)]);
    axis(gca, 'off'); % axis tight sets the axis limits to the range of the data.
    figtitleh=title([cmdData.fileName ' ' condition]);
    set(figtitleh,'Interpreter','none');
    
    %% Plot sdf
    sdfplot=subplot(2,1,2,'Layer','top','Parent', mainfig);
    %sdfh = axes('Position', [.15 .65 .2 .2], 'Layer','top');
    title('Spike Density Function','FontName','calibri','FontSize',11);
    hold on;
    if size(rasters(~isnantrial,:),1)<5 %if less than 5 good trials
        %useless plotting this
        sumall=NaN;
    else
        sumall=sum(rasters(~isnantrial,start-fsigma*3:stop+fsigma*3));
    end
    %     sdf=spike_density(sumall,fsigma)./length(find(~isnantrial)); %instead of number of trials
    sdf=fullgauss_filtconv(sumall,fsigma,causker)./length(find(~isnantrial)).*1000;
    %     sdf=sdf(fsigma+1:end-fsigma);
    
    %% calculate confidence intervals
    lcut_rasters=rasters(~isnantrial,start-fsigma*3:stop+fsigma*3);
    smoothtrial=zeros(size(lcut_rasters,1),size(lcut_rasters,2)-fsigma*6);
    for crsem=1:size(rasters(~isnantrial),1)
        smoothtrial(crsem,:)=fullgauss_filtconv(lcut_rasters(crsem,:),fsigma,causker).*1000;
    end
    %     smoothtrial=smoothtrial(:,fsigma+1:end-fsigma);
    if numrast==2 && rastnum==1  %collect first batch of trials
        first_smtrials=smoothtrial;
    end
    rastsem=std(smoothtrial)/ sqrt(size(smoothtrial,1)); %standard error of the mean
    %norminv([.025 .975], mean(smoothtrial), std(smoothtrial));
    rastsem = rastsem * 1.96; % 95% of the data will fall within 1.96 standard deviations of a normal distribution
    
    % testif significant diff
    %         differential spike density
    %  function exceeded by 2 SD the mean difference in activity
    %  during the 600-ms interval before target presentation, pro
    % vided the difference reached 6 SD and remained >2 SD
    % threshold for 50 ms.
    
    %     if numrast==2 && rastnum==numrast
    %         diff_trials=mean(first_smtrials)-mean(smoothtrial);
    %         diff_preal_epoch=diff_trials(1:alignidx-start); %alignidx-start-200
    %         difftime_preal=find(abs(diff_preal_epoch)>2*(std(diff_preal_epoch)),1);
    %         diff_postal_epoch=diff_trials(alignidx-start:alignidx-start+200);
    %         difftime_postal=find(abs(diff_postal_epoch)>2*(std(diff_postal_epoch)),1);
    %         if ~isempty(difftime_preal)
    %             %recursive time search
    %             difftime_preal=difftime_preal-find(abs(diff_trials(alignidx-start-200+difftime_preal+1:-1:1))<=2*(std(diff_preal_epoch)),1);
    %             difftime_preal=alignidx-start-200+1+difftime_preal;
    %         end
    %         if ~isempty(difftime_postal)
    %             %recursive time search
    %             difftime_postal=difftime_postal-find(abs(diff_trials(alignidx-start+difftime_postal+1:-1:1))<=2*(std(diff_preal_epoch)),1);
    %             difftime_postal=alignidx-start+1+difftime_postal;
    %         end
    %     else
    %         difftime_preal=[];
    %         difftime_postal=[];
    %     end
    if size(rasters(~isnantrial,:),1)>=5
        %    plot confidence intervals
        patch([1:length(sdf),fliplr(1:length(sdf))],[sdf-rastsem,fliplr(sdf+rastsem)],cmap(rastnum,:),'EdgeColor','none','FaceAlpha',0.1);
        %plot sdf
        plot(sdf,'Color',cmap(rastnum,:),'LineWidth',1.8);
        
        %         if ~isempty(difftime_preal)
        %             plot(difftime_preal,max([sdf(difftime_preal)-40 1]),'r*')
        %         end
        %         if ~isempty(difftime_postal)
        %             plot(difftime_postal,max([sdf(difftime_postal)-40 1]),'r*')
        %         end
    end
    
    % axis([0 stop-start 0 200])
    axis(gca,'tight');
    box off;
    set(gca,'Color','white','TickDir','out','FontName','calibri','FontSize',8); %'YAxisLocation','rigth'
    set(gca,'XTick',1:100:length(sdf),'XTickLabel',-plotstart:100:plotstop,'TickDir','out','box','off'); %
    
    %     hxlabel=xlabel(gca,'Time (ms)','FontName','calibri','FontSize',8);
    %     set(hxlabel,'Position',get(hxlabel,'Position') - [180 -0.2 0]); %doesn't stay there when export !
    hylabel=ylabel(gca,'Firing rate (spikes/s)','FontName','Calibri','FontSize',8);
    
    %     %% Plot eye velocities
    %     heyevelplot=subplot(2,1,(2*2/3)+1:2,'Layer','top','Parent', mainfig);
    %     title('Mean Eye Velocity','FontName','calibri','FontSize',11);
    %     hxlabel=xlabel(gca,'Time (ms)','FontName','calibri','FontSize',8);
    %
    %     hold on;
    %     if ~isempty(rasters)
    %         eyevel=alignedata(rastnum).eyevel;
    %         eyevel=mean(eyevel(:,start:stop));
    %         heyevelline(rastnum)=plot(eyevel,'Color',cmap(rastnum,:),'LineWidth',1.8);
    %         axis(gca,'tight');
    %         %         eyevelymax=max(eyevel);
    %         %         if eyevelymax>0.8
    %         %             eyevelymax=eyevelymax*1.1;
    %         %         else
    %         %             eyevelymax=0.8;
    %         %         end
    %         %         axis([0 stop-start 0 eyevelymax]);
    %         set(gca,'Color','none','TickDir','out','FontSize',8,'FontName','calibri','box','off');
    %         set(gca,'XTick',1:100:length(sdf),'XTickLabel',-plotstart:100:plotstop,'TickDir','out','box','off'); %
    %         ylabel(gca,'Eye velocity (deg/ms)','FontName','calibri','FontSize',8);
    %
    %         % get directions for the legend
    %         if isfield(alignedata,'dir')
    %             curdir{rastnum}=alignedata(rastnum).dir;
    %         else
    %             curdir{rastnum}='no';
    %             %             % need eyeh and eyev. But see cmdData.allprefdir
    %             %             sacdeg=nan(size(alignedata(1,rastnum).trialnb,2),1);
    %             %             for eyetr=1:size(alignedata(1,rastnum).trialnb,2)
    %             %                 thissach=alignedata(1,rastnum).eyeh(eyetr,alignedata(1,rastnum).alignt:alignedata(1,rastnum).alignt+100);
    %             %                 thissacv=alignedata(1,rastnum).eyev(eyetr,alignedata(1,rastnum).alignt:alignedata(1,rastnum).alignt+100);
    %             %                 minwidth=5;
    %             %                 [~, ~, thissacvel, ~, ~, ~] = cal_velacc(thissach,thissacv,minwidth);
    %             %                 peakvel=find(thissacvel==max(thissacvel),1);
    %             %                 sacendtime=peakvel+find(thissacvel(peakvel:end)<=...
    %             %                     (min(thissacvel(peakvel:end))+(max(thissacvel(peakvel:end))-min(thissacvel(peakvel:end)))/10),1);
    %             %                 try
    %             %                     sacdeg(eyetr)=abs(atand((thissach(sacendtime)-thissach(1))/(thissacv(sacendtime)-thissacv(1))));
    %             %                 catch
    %             %                     thissacv;
    %             %                 end
    %             %
    %             %                 % sign adjustements
    %             %                 if thissacv(sacendtime)<thissacv(1) % negative vertical amplitude -> vertical flip
    %             %                     sacdeg(eyetr)=180-sacdeg(eyetr);
    %             %                 end
    %             %                 if thissach(sacendtime)>thissach(1)%inverted signal: leftward is in postive range. Correcting to negative.
    %             %                     sacdeg(eyetr)=360-sacdeg(eyetr); % mirror image;
    %             %                 end
    %             %             end
    %             %             % a quick fix to be able to put "upwards" directions together
    %             %             distrib=hist(sacdeg,3); %floor(length(sacdeg)/2)
    %             %             if max(bwlabel(distrib,4))>1 && distrib(1)>1 && distrib(end)>1 %=bimodal distribution with more than 1 outlier
    %             %                 sacdeg=sacdeg+45;
    %             %                 sacdeg(sacdeg>360)=-(360-(sacdeg(sacdeg>360)-45));
    %             %                 sacdeg(sacdeg>0)= sacdeg(sacdeg>0)-45;
    %             %             end
    %             %             sacdeg=abs(median(sacdeg));
    %             %
    %             %             if sacdeg>45/2 && sacdeg <= 45+45/2
    %             %                 curdir{rastnum}='up_right';
    %             %             elseif sacdeg>45+45/2 && sacdeg <= 90+45/2
    %             %                 curdir{rastnum}='rightward';
    %             %             elseif sacdeg>90+45/2 && sacdeg <= 135+45/2
    %             %                 curdir{rastnum}='down_right';
    %             %             elseif sacdeg>135+45/2 && sacdeg < 180+45/2
    %             %                 curdir{rastnum}='downward';
    %             %             elseif sacdeg>=180+45/2 && sacdeg <= 225+45/2
    %             %                 curdir{rastnum}='down_left';
    %             %             elseif sacdeg>225+45/2 && sacdeg <= 270+45/2
    %             %                 curdir{rastnum}='leftward';
    %             %             elseif sacdeg>270+45/2 && sacdeg <= 315+45/2
    %             %                 curdir{rastnum}='up_left';
    %             %             else
    %             %                 curdir{rastnum}='upward';
    %             %             end
    %         end
    %         aligntype{rastnum}=alignedata(rastnum).alignlabel;
    %     else
    %         curdir{rastnum}='no';
    %         aligntype{rastnum}='data';
    %     end
end

if ~isempty(rasters)
    % drawing the alignment bars
    axes(sdfplot)
    currylim=get(gca,'YLim');
    patch([repmat((alignidx-start)-2,1,2) repmat((alignidx-start)+2,1,2)], ...
        [[0 currylim(2)] fliplr([0 currylim(2)])], ...
        [0 0 0 0],[0 0 1],'EdgeColor','none','FaceAlpha',0.5);
    %     axes(heyevelplot)
    %     currylim=get(gca,'YLim');
    %     patch([repmat((alignidx-start)-2,1,2) repmat((alignidx-start)+2,1,2)], ...
    %         [[0 currylim(2)] fliplr([0 currylim(2)])], ...
    %         [0 0 0 0],[0 0 1],'EdgeColor','none','FaceAlpha',0.5);
end

if singleSSD
    % plot SSD bar
    axes(sdfplot)
    currylim=get(gca,'YLim');
    alignTime=alignidx-start+prevssd;
    patch([repmat((alignTime)-2,1,2) repmat((alignTime)+2,1,2)], ...
        [[0 currylim(2)] fliplr([0 currylim(2)])], ...
        [0 0 0 0],[1 0 0],'EdgeColor','none','FaceAlpha',0.5);
    % plot SSRT
    alignTime=alignTime+cmdData.allmssrt_tacho{1, 1}{1};
    patch([repmat((alignTime)-2,1,2) repmat((alignTime)+2,1,2)], ...
        [[0 currylim(2)] fliplr([0 currylim(2)])], ...
        [0 0 0 0],[.5 .5 .5],'EdgeColor','none','FaceAlpha',0.5);
end

% Save figure

figDir='E:\Dropbox\Vincent Docs\CbTimingPredict\figures\SingleNeuronExample\';
exportfigname=[figDir cmdData.fileName 'targetalign'];
%     print(exportfig, '-dpng', '-noui', '-opengl','-r600', exportfigname);
print(gcf, '-dsvg', '-noui', '-painters','-r600', exportfigname);
delete(gcf);

