function CmdBehaviorResults(gsdata,conn)
%plot saccade latencies, psychometric and tachometric curves
global directory slash;
if isempty(directory)
    [directory,slash]=SetUserDir;
end

unitsDBinfo=gsdata.alldb;
unitList=cellfun(@(x) x.unit_id, unitsDBinfo);
% get cluster index from database
dbConn = connect2DB('vp_sldata');
clusterIdx=getclusterindex(dbConn,unitList);

%number cells
gs.goodrecs=~(isnan(clusterIdx) | cellfun('isempty',gsdata.allsacdelay));

%remove bad apples
fn = fieldnames(gsdata);
for lp=1:length(fn)
    gsdata.(fn{lp})=gsdata.(fn{lp})(gs.goodrecs,:);
end

gs.cellnum=sum(~cellfun('isempty',gsdata.allsacdelay));

disp([num2str(gs.cellnum) ' cells for countermanding task ' ])

%% individual sessions

%% saccade delays
allsacdelays=cellfun(@(x) x.nsst,gsdata.allsacdelay,'UniformOutput',false);
allsacdelays=[allsacdelays{:}];
xlim=ceil(max(allsacdelays)/100)*100;
[saclatquant,saclatxlims]=hist(allsacdelays,0:25:xlim);
saclatfreq=saclatquant./sum(saclatquant);
%
% figure;plot(saclatxlims,gauss_filtconv(saclatfreq,1));
% % title('Saccade Latency','FontName','calibri','FontSize',15);
% hxlabel=xlabel(gca,'Saccade latency','FontName','calibri','FontSize',12);
% set(gca,'XTick',0:200:max(get(gca,'xlim')),'TickDir','out','box','off');
% hylabel=ylabel(gca,'Proportion of saccades','FontName','calibri','FontSize',12);
% curylim=get(gca,'YLim');
% set(gca,'TickDir','out','box','off');

%% psychometric and tachometric curves

query = ['SELECT a_file FROM sorts s INNER JOIN recordings r on s.recording_id_fk = r.recording_id WHERE sort_id IN (' ...
    sprintf('%.0f,' ,cellfun(@(x) x.sort_id,gsdata.alldb(1:end-1,1))) num2str(gsdata.alldb{end,1}.sort_id) ')'];
gs.recnames=fetch(conn,query);

behavData=struct('subject',[],'mssrt',[],'inhibfun',[],'ccssd',[],'nccssd',[],...
    'ssds',[],'prevssds',[],'trialidx',[],'alltacho',[]);

for gsFileNum=1:size(gs.recnames,1)
    
    gs.recnames{gsFileNum,1}=gs.recnames{gsFileNum,1}(1:end-1);
    
    %% subject and procdir
    if strcmp('R',gs.recnames{gsFileNum,1}(1))
        behavData(gsFileNum).subject='Rigel';
        procdir = [directory,'processed',slash,'Rigel',slash];
    elseif strcmp('S',gs.recnames{gsFileNum,1}(1))
        behavData(gsFileNum).subject='Sixx';
        procdir = [directory,'processed',slash,'Sixx',slash];
    elseif strcmp('H',gs.recnames{gsFileNum,1}(1))
        behavData(gsFileNum).subject='Hilda';
        procdir = [directory,'processed',slash,'Hilda',slash];
    end
    
    %% what file to load
    procdirlisting=dir(procdir);
    procdirfileNames={procdirlisting.name};
    loadfile=procdirfileNames(~cellfun('isempty',regexpi(procdirfileNames,gs.recnames{gsFileNum,1},'match')));
    
    % if two versions (REX and Spike2), chose REX (for now)
    if length(loadfile)>1
        loadfile=loadfile(~cellfun('isempty',regexpi(loadfile,'REX','match')));
        if isempty(loadfile) %actually it's an old file that sneaked in
            loadfile=procdirfileNames(~cellfun('isempty',regexpi(procdirfileNames,gs.recnames{gsFileNum,1},'match')));
            loadfile=loadfile(~cellfun('isempty',regexpi(loadfile,'REX','match')));
        end
    end
    [mssrt,inhibfun,ccssd,nccssd,ssds,...
        ~,~,~,~,prevssds,trialidx,alltacho]=findssrt(loadfile{:});
    behavData(gsFileNum).mssrt=mssrt;
    behavData(gsFileNum).inhibfun=inhibfun;
    behavData(gsFileNum).ccssd=ccssd;
    behavData(gsFileNum).nccssd=nccssd;
    behavData(gsFileNum).ssds=ssds;
    behavData(gsFileNum).prevssds=prevssds;
    behavData(gsFileNum).trialidx=trialidx;
    behavData(gsFileNum).alltacho=alltacho;
    clear mssrt inhibfun ccssd nccssd ssds prevssds trialidx alltacho
end

xtachLims=cellfun(@(x) [nanmin(x.xtach) nanmax(x.xtach)], {behavData.alltacho},'UniformOutput',false);
xtachComp=nan(size(gs.recnames,1),nanmax([xtachLims{:}])-nanmin([xtachLims{:}])+1); 
tachofig=figure('color','white','position',[1017 184 560 420]); cmap=colormap('lines');
tachoplots=subplot(1,2,1); hold on;
for gsFileNum=1:size(gs.recnames,1)
    if ~isempty(behavData(gsFileNum).alltacho.tachowidth) && behavData(gsFileNum).alltacho.tachowidth>0
        xtach=behavData(gsFileNum).alltacho.xtach;
        tach=behavData(gsFileNum).alltacho.tach;
        filttach=gauss_filtconv(tach,4);
        coretach=filttach(xtach>-20 & xtach<170);
        %         coretach(isnan(coretach))=0;
        plot(tachoplots,xtach(xtach>-20 & xtach<170),coretach);
        xtachComp(gsFileNum,xtach(1)-min([xtachLims{:}])+1:xtach(end)-min([xtachLims{:}])+1)=filttach;
    end
end

%% Tachometric curves
title('Tachometric curve','FontName','calibri','FontSize',12);
hxlabel=xlabel(gca,'rPT (ms)','FontName','calibri','FontSize',12);
set(gca,'Xlim',[-20 170],'XTick',-20:50:170,'TickDir','out','box','off'); %'XTickLabel',[50:50:400]
hylabel=ylabel(gca,'Fraction cancelled','FontName','calibri','FontSize',12);
set(gca,'Ylim',[0 1],'TickDir','out','box','off');

medtachoplot=subplot(1,2,2);
meanTach=nanmean(xtachComp);
meanTach=meanTach(-19-min([xtachLims{:}]): 171-min([xtachLims{:}]));
meanTachsem=(nanstd(xtachComp)/ sqrt(size(xtachComp,1))) * 1.96;
meanTachsem=meanTachsem(-19-min([xtachLims{:}]): 171-min([xtachLims{:}]));
% medTachsem(isnan(medTachsem))=0;
plot(meanTach,'linewidth',2)
patch([1:length(meanTachsem),fliplr(1:length(meanTachsem))],...
    [meanTach-meanTachsem,fliplr(meanTach+meanTachsem)],cmap(1,:),'EdgeColor','none','FaceAlpha',0.5);
% patch([1:11,fliplr(1:11)],[1:11,fliplr(5:15)],cmap(1,:),'EdgeColor','none','FaceAlpha',0.5);

title('Mean tachometric curve','FontName','calibri','FontSize',12);
xlabel(gca,'rPT (ms)','FontName','calibri','FontSize',12);
set(gca,'Xlim',[0 190],'XTick',0:50:190,'XTickLabel',-20:50:170,'TickDir','out','box','off');
ylabel(gca,'Fraction cancelled','FontName','calibri','FontSize',12);
set(gca,'Ylim',[0 1],'TickDir','out','box','off');

% plot for each individual
subjectList={behavData.subject};
subjects={'Rigel';'Sixx';'Hilda'};
meanTach=nanmean(xtachComp);
meanTach=meanTach(-19-min([xtachLims{:}]): 171-min([xtachLims{:}]));
meanTachsem=(nanstd(xtachComp)/ sqrt(size(xtachComp,1))) * 1.96;
meanTachsem=meanTachsem(-19-min([xtachLims{:}]): 171-min([xtachLims{:}]));

figure; hold on;
plot(meanTach,'linewidth',2)
patch([1:length(meanTachsem),fliplr(1:length(meanTachsem))],...
    [meanTach-meanTachsem,fliplr(meanTach+meanTachsem)],cmap(1,:),'EdgeColor','none','FaceAlpha',0.5);

for subjNum=1:3
    indivxtachComp=xtachComp(cellfun(@(x) strcmp(x,subjects{subjNum}), subjectList),:);
    meanTach=nanmean(indivxtachComp);
    meanTach=meanTach(-19-min([xtachLims{:}]): 171-min([xtachLims{:}]));
    meanTachsem=(nanstd(indivxtachComp)/ sqrt(size(indivxtachComp,1))) * 1.96;
    meanTachsem=meanTachsem(-19-min([xtachLims{:}]): 171-min([xtachLims{:}]));
    plot(meanTach,'linewidth',2)
    if sum(meanTachsem)>0
        patch([1:length(meanTachsem),fliplr(1:length(meanTachsem))],...
            [meanTach-meanTachsem,fliplr(meanTach+meanTachsem)],cmap(subjNum+1,:),'EdgeColor','none','FaceAlpha',0.5);
    end
end
legend({'All','','Rigel','Sixx','','Hilda',''})

title('Mean tachometric curve per individual','FontName','calibri','FontSize',12);
xlabel(gca,'rPT (ms)','FontName','calibri','FontSize',12);
set(gca,'Xlim',[0 190],'XTick',0:50:190,'XTickLabel',-20:50:170,'TickDir','out','box','off');
ylabel(gca,'Fraction cancelled','FontName','calibri','FontSize',12);
set(gca,'Ylim',[0 1],'TickDir','out','box','off');

%% rPT figures
[rPTcComp,rPTeComp]=deal(nan(size(gs.recnames,1),max([xtachLims{:}])-min([xtachLims{:}])+1));
rPTfig=figure('color','white','position',[1054 -170 560 420]);
rPTplots=subplot(1,2,1); hold on;
for gsFileNum=1:size(gs.recnames,1)
    if ~isempty(behavData(gsFileNum).alltacho.tachowidth) && behavData(gsFileNum).alltacho.tachowidth>0
        xtach=behavData(gsFileNum).alltacho.xtach;
        rPTc=behavData(gsFileNum).alltacho.rPTc;
        rPTe=behavData(gsFileNum).alltacho.rPTe;
        rPTce=rPTc+rPTe;
        filtrPTc = gauss_filtconv(rPTc./max(rPTce),10);
        filtrPTe = gauss_filtconv(rPTe./max(rPTce),10);
        plot(xtach,filtrPTc,'color',cmap(1,:));
%         hold on
        %           plot(xtach,rPTc,'r')
        %           plot(xtach,rPTe,'b')
        plot(xtach,filtrPTe,'color',cmap(2,:));
        rPTcComp(gsFileNum,xtach(1)-min([xtachLims{:}])+1:xtach(end)-min([xtachLims{:}])+1)=filtrPTc;
        rPTeComp(gsFileNum,xtach(1)-min([xtachLims{:}])+1:xtach(end)-min([xtachLims{:}])+1)=filtrPTe;
    end
end
title('Distribution of rPTs','FontName','calibri','FontSize',15);
hxlabel=xlabel(gca,'rPT (ms)','FontName','calibri','FontSize',12);
set(gca,'Xlim',[-100 400],'XTick',-100:50:400,'XTickLabel',-100:50:400); 
hylabel=ylabel(gca,'Normalized frequency','FontName','calibri','FontSize',12);
set(gca,'TickDir','out','box','off','ylim',[0 1]);
% set(gca,'XTick',-200:100:500);
% set(gca,'XTickLabel',-200:100:500);

meanrPTplot=subplot(1,2,2); hold on;
meanrPTc=nanmean(rPTcComp); meanrPTe=nanmean(rPTeComp);
meanrPTce=meanrPTc+meanrPTe;
plot(meanrPTc(-99-min([xtachLims{:}]): 401-min([xtachLims{:}]))./max(meanrPTce),'linewidth',2)
plot(meanrPTe(-99-min([xtachLims{:}]): 401-min([xtachLims{:}]))./max(meanrPTce),'linewidth',2)
set(gca,'Xlim',[0 500],'XTick',0:50:500,'XTickLabel',-100:50:400); 
title('Distribution of mean rPTs','FontName','calibri','FontSize',15);
xlabel(gca,'rPT (ms)','FontName','calibri','FontSize',12);
ylabel(gca,'Normalized frequency','FontName','calibri','FontSize',12);
set(gca,'TickDir','out','box','off','ylim',[0 1]);

%% SSRTs
SSRTs=[behavData.mssrt];
figure;
subjects={'Rigel';'Sixx';'Hilda'};
for subjNum=1:3
%     subjects{subjNum}
    indivSSRTs=SSRTs(cellfun(@(x) strcmp(x,subjects{subjNum}), {behavData.subject}));
    indivSSRTs=indivSSRTs(~isnan(indivSSRTs));
    subplot(3,1,subjNum);
    histogram(indivSSRTs,50:10:130);
    title(['SSRTs for ' subjects{subjNum}])
%     numel(indivSSRTs)
%     meanIndivSSRT=mean(indivSSRTs)
%     SEIndivSSRT=std(indivSSRTs)/sqrt(numel(indivSSRTs))
end
figure;
histogram(SSRTs,50:10:130);
title('SSRTs, all subjects')
% numel(SSRTs)
% meanAllSSRT=mean(SSRTs)
% SEAllSSRT=std(SSRTs)/sqrt(numel(SSRTs))

% save('behavData','behavData');

%% psychometric curves
allinhibFuns={behavData(~cellfun(@(x) isnan(sum(x)),{behavData.inhibfun})).inhibfun};
allSSDs={behavData(~cellfun(@(x) isnan(sum(x)),{behavData.ssds})).ssds};
% xInfunLims=cellfun(@(x) [nanmin(x) nanmax(x)], allSSDs,'UniformOutput',false);
% max(cellfun(@(x) x(end)-x(1), allSSDs))
inhibFunsComp=nan(size(allSSDs,2),20); 
inhibfunfig=figure('color','white','position',[1017 184 560 420]); cmap=colormap('lines');
inhibfunplots=subplot(1,2,1); hold on;
for gsFileNum=1:size(allSSDs,2)
        SSDS=allSSDs{gsFileNum};
        inhibFun=allinhibFuns{gsFileNum};
        plot(SSDS,inhibFun);
    %         change into quartiles !
    % interpolate
    %add time points to make inhibition function spread over 250ms, from 0 to 1
    timeoutliers=(250-(SSDS(end)-SSDS(1)))/2;
    rs_inhibFun = timeseries([0;inhibFun;1],[SSDS(1)-timeoutliers; SSDS; SSDS(end)+timeoutliers]);
    rs_inhibFun.TimeInfo.Units='milliseconds';
    rs_inhibFun = resample(rs_inhibFun, linspace(SSDS(1)-timeoutliers,SSDS(end)+timeoutliers,20));
    % figure; plot(ts1)
    inhibFunsComp(gsFileNum,:)=rs_inhibFun.Data;
end

subplot(1,2,2);
meanInhibFun=nanmean(inhibFunsComp);meanInhibFun=fliplr(meanInhibFun);
meanInhibFunsem=(nanstd(inhibFunsComp)/ sqrt(size(inhibFunsComp,1))) * 1.96;
meanInhibFunsem=fliplr(meanInhibFunsem);
plot(meanInhibFun,'linewidth',2)
patch([1:length(meanInhibFunsem),fliplr(1:length(meanInhibFunsem))],...
    [meanInhibFun-meanInhibFunsem,fliplr(meanInhibFun+meanInhibFunsem)],cmap(1,:),'EdgeColor','none','FaceAlpha',0.5);
% patch([1:11,fliplr(1:11)],[1:11,fliplr(5:15)],cmap(1,:),'EdgeColor','none','FaceAlpha',0.5);

title('Mean inhibtion function', 'FontName','calibri','FontSize',12);
xlabel(gca,'SSDS(ms)','FontName','calibri','FontSize',12);
set(gca,'Xlim',[0 20],'XTick',0:5:20,'XTickLabel',round(linspace(-125,125,5)),'TickDir','out','box','off');
ylabel(gca,'Fraction cancelled','FontName','calibri','FontSize',12);
set(gca,'Ylim',[0 1],'TickDir','out','box','off');

% plot for each individual
subjectList={behavData(~cellfun(@(x) isnan(sum(x)),{behavData.inhibfun})).subject};
subjects={'Rigel';'Sixx';'Hilda'};

figure; hold on;
plot(meanInhibFun,'linewidth',2)
patch([1:length(meanInhibFunsem),fliplr(1:length(meanInhibFunsem))],...
    [meanInhibFun-meanInhibFunsem,fliplr(meanInhibFun+meanInhibFunsem)],cmap(1,:),'EdgeColor','none','FaceAlpha',0.5);

for subjNum=1:3
    indivInhibFunsComp=inhibFunsComp(cellfun(@(x) strcmp(x,subjects{subjNum}), subjectList),:);
    meanInhibFun=nanmean(indivInhibFunsComp);meanInhibFun=fliplr(meanInhibFun);
    meanInhibFunsem=(nanstd(indivInhibFunsComp)/ sqrt(size(indivInhibFunsComp,1))) * 1.96;
    plot(meanInhibFun,'linewidth',2,'color',cmap(subjNum+1,:));
    patch([1:length(meanInhibFunsem),fliplr(1:length(meanInhibFunsem))],...
        [meanInhibFun-meanInhibFunsem,fliplr(meanInhibFun+meanInhibFunsem)],cmap(subjNum+1,:),'EdgeColor','none','FaceAlpha',0.5);
end
legend({'All','','Rigel','','Sixx','','Hilda',''})
title('Mean inhibition function per individual', 'FontName','calibri','FontSize',12);
xlabel(gca,'SSDS(ms)','FontName','calibri','FontSize',12);
set(gca,'Xlim',[0 20],'XTick',0:5:20,'XTickLabel',round(linspace(-125,125,5)),'TickDir','out','box','off');
ylabel(gca,'Fraction cancelled','FontName','calibri','FontSize',12);
set(gca,'Ylim',[0 1],'TickDir','out','box','off');


