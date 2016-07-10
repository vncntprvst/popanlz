function sing_sess_countermanding(data)
% Analyze and display data from a single session of countermanding task

global directory slash;

%% settings
userinfo=SetUserDir;
directory=userinfo.directory;
slash=userinfo.slash;

if isstr(data)
%% finding best candidate

% userinfo.user,userinfo.dbldir,userinfo.mapdr,userinfo.servrep,userinfo.mapddataf
CCNdb = connect2DB('vp_sldata');

cd(userinfo.syncdir);
load('cDn_gsdata.mat'); %cDn_gsdata.mat  top_cortex_gsdata.mat

%number cells
goodrecs=~cellfun('isempty',gsdata.allsacdelay);
% st.goodrecs=~cellfun('isempty',stdata.allsacdelay);

%remove bad apples
fn = fieldnames(gsdata);
for lp=1:length(fn)
    gsdata.(fn{lp})=gsdata.(fn{lp})(goodrecs,:);
end

sigma=15;
periEventActivity=cellfun(@(x) conv_raster(x(1,1).rast,sigma,x(1,1).alignt-(300+sigma*3),x(1,1).alignt+(299+sigma*3)), gsdata.allndata(:,1), 'UniformOutput',false); %600ms epoch

% remove recordings with one trial
properLengthRecs=cellfun(@(x) length(x)>1,periEventActivity);
goodrecs(goodrecs)=properLengthRecs;
periEventActivity=periEventActivity(properLengthRecs);

% get unit cluster info and profiles
unit_ids=cellfun(@(x) x.unit_id,gsdata.alldb);

[sorted_unit_ids,sunitid_idx]=sort(unit_ids);
query = ['SELECT c.profile, c.profile_type FROM clusters c WHERE cluster_id IN (' sprintf('%.0f,' ,sorted_unit_ids(1:end-1)') num2str(sorted_unit_ids(end)) ')'];
profiles = fetch(CCNdb,query);
sunitid_revidx(sunitid_idx)=1:length(unit_ids);
clusidx=[profiles{sunitid_revidx,2}];
clustypes={profiles{sunitid_revidx,1}};

% find most active cells, aligned to saccade, in cluster #1
periEventActivity=periEventActivity(clusidx==101);
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
numTrials=cellfun(@(x) [size(x(1).rast,1) size(x(2).rast,1) size(x(3).rast,1)],gsdata.allndata(:,2), 'UniformOutput',false);
% mostTrials=cat(1,numTrials{[14,5,4,31,15]});
mostTrials=cat(1,numTrials{:});
mostTrials=mostTrials(clusidx==101,:);
[~,mostNSST]=sort(mostTrials(:,1),'descend'); [~,mostNSSTscore]=sort(mostNSST);
[~,mostCSST]=sort(mostTrials(:,2),'descend'); [~,mostCSSTscore]=sort(mostCSST);
[~,mostNCSST]=sort(mostTrials(:,3),'descend'); [~,mostNCSSTscore]=sort(mostNCSST);
[~,mostTrials]=sort(sum([mostNSSTscore,mostCSSTscore,mostNCSSTscore],2));[~,mostTrials]=sort(mostTrials);

%combine two indices to find recordings with best activity and most trials
[~,bestFiles]=sort(topEventActivity+mostTrials);

% hand picked recordings with best activity and most trials were 5,14,15

% plot them
% for cellNum=1:10
%     bestFile=bestFiles(cellNum);
%     
%     figure('position',[1686 49 867 1308])
%     for plotNum=1:3
%         %get rasters
%         periTargetActivity=cellfun(@(x) conv_raster(x(plotNum).rast,sigma,x(plotNum).alignt-(100+sigma*3),x(plotNum).alignt+(500+sigma*3)), gsdata.allndata(:,2), 'UniformOutput',false); %600ms epoch
%         periTargetActivity=periTargetActivity(clusidx==101);
%         
%         subplot(3,1,plotNum)
%         plot(periTargetActivity{bestFile},'LineWidth',2);
%         currylim=get(gca,'YLim');
%             patch([repmat(98,1,2) repmat(102,1,2)], ...
%             [[0 currylim(2)] fliplr([0 currylim(2)])], ...
%             [0 0 0 0],[1 0 0],'EdgeColor','none','FaceAlpha',0.5);
%     
%     axis(gca,'tight');
%     box off;
%         
%     end
% end

%best cell is bestFiles(5) = 31 (in cluster #1)

end

%% session psychophysics

% subjects={'Rigel','Sixx','Hilda'};
% 
% % best use a recording made with gapstop training and lots of trials
% subject=subjects{3};
% recname='H48L5A1_18102'; % S148cnttrain
% load(recname,'allcodes','alltimes','allbad','saccadeInfo');
% % load([recname,'_sac'],'dataaligned');
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


%% session recording

    rasters=alignedata(rastnum).rasters;
    alignidx=alignedata(rastnum).alignidx;
    if chrono
       cut_chrasters=zeros(length([alignedata.trials]),plotstart+plotstop+1);
       %listing relevant trials in a continuous series with other rasters
       chronoidx=ismember(sort([alignedata.trials]),alignedata(rastnum).trials); 
    end
    greyareas=alignedata(rastnum).allgreyareas;
    start=alignidx - plotstart;
    stop=alignidx + plotstop;
    
    if start < 1
        start = 1;
    end
    if stop > length(rasters)
        stop = length(rasters);
    end
    
    %trials = size(rasters,1);
    isnantrial=zeros(1,size(rasters,1));
    
    if chrono
        onerastplot=subplot(numsubplot,1,1:(numsubplot/3),'Layer','top', ...
            'XTick',[],'YTick',[],'XColor','white','YColor','white', 'Parent', handles.mainfig);
    else
        if numrast==1
            hrastplot(rastnum)=subplot(numsubplot,1,1:2,'Layer','top', ...
                'XTick',[],'YTick',[],'XColor','white','YColor','white', 'Parent', handles.mainfig);
        else
            hrastplot(rastnum)=subplot(numsubplot,1,rastnum,'Layer','top', ...
                'XTick',[],'YTick',[],'XColor','white','YColor','white', 'Parent', handles.mainfig);
        end
    end
    %reducing spacing between rasters
    if numrast>1 && ~chrono
        rastpos=get(gca,'position');
        rastpos(2)=rastpos(2)+rastpos(4)*0.5;
        set(gca,'position',rastpos);
    end
    
    % sorting rasters according greytime
    viscuetimes=nan(size(greyareas,2),2);
    for grst=1:size(greyareas,2)
        viscuetimes(grst,:)=greyareas{grst}(1,:);
    end
    
    if ~chrono
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
        
        
            if(size(rasters,1) == 1)
                plot([indx;indx],[indy;indy+1],'color',cc(rastnum,:),'LineStyle','-'); % plot rasters
            else
                plot([indx';indx'],[indy';indy'+1],'color',cc(rastnum,:),'LineStyle','-'); % plot rasters
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
            patch(grxlims, grylims, [0 0 0], 'EdgeColor', 'none','FaceAlpha', 0.2)
        end

    set(gca,'xlim',[1 length(start:stop)]);
    axis(gca, 'off'); % axis tight sets the axis limits to the range of the data.
    
    
    %% Plot sdf
    sdfplot=subplot(numsubplot,1,(numsubplot/3)+1:(numsubplot/3)+(numsubplot/3),'Layer','top','Parent', handles.mainfig);
    %sdfh = axes('Position', [.15 .65 .2 .2], 'Layer','top');
    title('Spike Density Function','FontName','calibri','FontSize',11);
    hold on;
    if size(rasters(~isnantrial,:),1)<5 %if less than 5 good trials
        %useless plotting this
        sumall=NaN;
    else
        sumall=sum(rasters(~isnantrial,start:stop));
    end
    %     sdf=spike_density(sumall,fsigma)./length(find(~isnantrial)); %instead of number of trials
    sdf=fullgauss_filtconv(sumall,fsigma,causker)./length(find(~isnantrial)).*1000;
%     sdf=sdf(fsigma+1:end-fsigma);
    
    %% calculate confidence intervals
    lcut_rasters=rasters(~isnantrial,start:stop);
    smoothtrial=zeros(size(lcut_rasters));
    for crsem=1:size(rasters(~isnantrial),1)
        smoothtrial(crsem,:)=fullgauss_filtconv(lcut_rasters(crsem,:),fsigma,causker).*1000; 
    end
%     smoothtrial=smoothtrial(:,fsigma+1:end-fsigma);
    if numrast==2 && rastnum==1  %collect old trials
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

        if numrast==2 && rastnum==numrast
        diff_trials=mean(first_smtrials)-mean(smoothtrial);
        diff_preal_epoch=diff_trials(alignidx-start-200:alignidx-start);
        difftime_preal=find(abs(diff_preal_epoch)>2*(std(diff_preal_epoch)),1);
        diff_postal_epoch=diff_trials(alignidx-start:alignidx-start+200);
        difftime_postal=find(abs(diff_postal_epoch)>2*(std(diff_postal_epoch)),1);
            if ~isempty(difftime_preal)
                %recursive time search
                difftime_preal=difftime_preal-find(abs(diff_trials(alignidx-start-200+difftime_preal+1:-1:1))<=2*(std(diff_preal_epoch)),1);
                difftime_preal=alignidx-start-200+1+difftime_preal;
            end
            if ~isempty(difftime_postal)
                %recursive time search
                difftime_postal=difftime_postal-find(abs(diff_trials(alignidx-start+difftime_postal+1:-1:1))<=2*(std(diff_preal_epoch)),1);
                difftime_postal=alignidx-start+1+difftime_postal;
            end
        else
            difftime_preal=[];
            difftime_postal=[];
        end
  if size(rasters(~isnantrial,:),1)>=5      
     %    plot confidence intervals
    patch([1:length(sdf),fliplr(1:length(sdf))],[sdf-rastsem,fliplr(sdf+rastsem)],cc(rastnum,:),'EdgeColor','none','FaceAlpha',0.1);
    %plot sdf
    plot(sdf,'Color',cc(rastnum,:),'LineWidth',1.8);
    
     if ~isempty(difftime_preal)
         plot(difftime_preal,max([sdf(difftime_preal)-40 1]),'r*')
     end
     if ~isempty(difftime_postal)
         plot(difftime_postal,max([sdf(difftime_postal)-40 1]),'r*')
     end
  end
  
    % axis([0 stop-start 0 200])
    axis(gca,'tight');
    box off;
    set(gca,'Color','white','TickDir','out','FontName','calibri','FontSize',8); %'YAxisLocation','rigth'
    %     hxlabel=xlabel(gca,'Time (ms)','FontName','calibri','FontSize',8);
    %     set(hxlabel,'Position',get(hxlabel,'Position') - [180 -0.2 0]); %doesn't stay there when export !
    hylabel=ylabel(gca,'Firing rate (spikes/s)','FontName','calibri','FontSize',8);
    currylim=get(gca,'YLim');
    
    if ~isempty(rasters)
        % drawing the alignment bar
        patch([repmat((alignidx-start)-2,1,2) repmat((alignidx-start)+2,1,2)], ...
            [[0 currylim(2)] fliplr([0 currylim(2)])], ...
            [0 0 0 0],[1 0 0],'EdgeColor','none','FaceAlpha',0.5);
    end

