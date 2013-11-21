%% directory
global directory slash;
if strcmp(getenv('username'),'SommerVD') || strcmp(getenv('username'),'vp35')
    directory = 'C:\Data\Recordings\';
elseif  strcmp(getenv('username'),'DangerZone')
    directory = 'E:\Data\Recordings\';
else
    directory = 'B:\data\Recordings\';
end
slash = '\';
%% files
%dentate
% CmdFileName={'S113L4A5_13500';'S114L4A5_14321';'S122L4A5_14010';'S126A4L6_13690';...
%     'R132L4P4_20152';'S121L4A5_13091';'S112l4a5_12971';...
%     'S117L4A6_12741';'S118L4A5_13081';'S115L4A6_12871';'S116L4A6_15431';...
%     'H62L5A5_22211';'S96L2A5_13651';'S99L2A5_13242';'H56L5A5_21502';...
%     'S116L4A6_15450';'H25L5A5_22760';'H53L5A5_20901';'S123L4A6_13691';...
%     'S125L4A6_13990';'R167L1A3_20671';'R166L1A3_20371';'S120L4A5_12912';...
%     'S105L1A8_16050';...
%     'S101L2A5_11251';'S89L4A5_12401';'S102L3A4_13001';'S111L4A5_13502';...
%     'H56L5A5_21321';'R97L7A1_19001'};

%Pause rebound cells
% CmdFileName={'S119L4A5_14391';'S111L4A5_14392';'S114L4A5_13650';...
%     'S121L4A5_13091';'R132L4P4_20152'}; % 'S105L1A8_16050' %where is that? ;'R97L7A1_19001';'S89L4A5_12401' %too different

% % %Ramping cells
CmdFileName={'S125L4A6_13990';'S123L4A6_13691';'H53L5A5_20901';...
    'H56L5A5_21502';'S96L2A5_13651';'S115L4A6_12871';...
    'S118L4A5_13081'}; % 'S114L4A5_14321' FEF like % 'S116L4A6_15450' 'S120L4A5_12912' 'R167L1A3_20671'

% % % poster boyz   
% CmdFileName={'S119L4A5_14391';'S123L4A6_13691'};

%For cancellation, for the moment keep 
% CmdFileName={'S96L2A5_13651';'S115L4A6_12871'};
%+ try to fix 'H53L5A5_20901'

%For conflict 'S114L4A5_13650';
%  CmdFileName={'S119L4A5_14391'};

%Fastigial & Vermis
%CmdFileName={'R167L1A3_20671';'R166L1A3_20371'}

for FileNb=1:length(CmdFileName);
    try
%% subject and procdir
if strcmp('R',CmdFileName{FileNb}(1))
   subject='Rigel';
   procdir = [directory,'processed',slash,'Rigel',slash];
elseif strcmp('S',CmdFileName{FileNb}(1))
   subject='Sixx';
   procdir = [directory,'processed',slash,'Sixx',slash];
elseif strcmp('H',CmdFileName{FileNb}(1))
   subject='Hilda';
   procdir = [directory,'processed',slash,'Hilda',slash];
end

procdirlisting=dir(procdir);
procdirfileNames={procdirlisting.name};
loadfile=procdirfileNames(~cellfun('isempty',regexpi(procdirfileNames,CmdFileName{FileNb},'match')));
% if two versions (REX and Spike2), chose Spike2
if length(loadfile)>1
    loadfile=loadfile(~cellfun('isempty',regexpi(loadfile,'Sp2','match')));
    if isempty(loadfile) %actually it's an old file that sneaked in
        loadfile=procdirfileNames(~cellfun('isempty',regexpi(procdirfileNames,CmdFileName{FileNb},'match')));
        loadfile=loadfile(~cellfun('isempty',regexpi(loadfile,'REX','match')));
    end
end

%% get countermanding session results 
[mssrt,inhibfun,ccssd,nccssd,ssdvalues,tachomc,tachowidth,sacdelay,rewtimes]=findssrt(loadfile{:}, 0);
mssrt=max([mssrt tachomc+tachowidth/2]);

%% align rasters
% set presets
tasktype='gapstop';
[~, trialdirs] = data_info(loadfile{:}, 1, 1); %reload file: yes (shouldn't happen, though), skip unprocessed files: yes

%alignments=1:3;

% sac vs stop
% firstalign=6; 
% secondalign=8; 
% aligntype='failed_fast';
% triplot = 0; 
%     plotstart=1000;
%     plotstop=1000;
%  [ctmatchlatidx nctmatchlatidx]=deal(NaN);
    
% tgt vs stop
firstalign=7; 
secondalign=8; 
aligntype='correct_slow';
plottype = 3; % 3 for splitting data in three groups: short SSD, med SSD and long SSD
plotstart=200;
plotstop=600;
 [ctmatchlatidx nctmatchlatidx]=deal(NaN);

% ssd
% firstalign=7; % as if align to target
% secondalign=507; 
% aligntype='ssd';
% triplot = 0; 
%     plotstart=800;
%     plotstop=600;

includebad=0;
spikechannel=1; %select appropriate cluster 
keepdir='compall'; %alldir %which sac directions
togrey=[];
singlerastplot=0;

if strcmp(aligntype,'ssd')
% if aligning to ssd, got to align NSS trials according to latency 
%allssds=unique([ccssd;nccssd]);

% canceled trials
ccssdval=unique(ccssd);
ctmatchlatidx=zeros(length(sacdelay),length(ccssdval));
for ssdval=1:length(ccssdval)
ctmatchlatidx(:,ssdval)=sacdelay>ccssdval(ssdval)+round(mssrt);
end
nullidx=sum(ctmatchlatidx,2)==0;
ctmatchlatidx(nullidx,1)=1;
% getting ssds for each NNS trial, taking the highest ssd. 
ctmatchlatidx=ccssdval(sum(ctmatchlatidx,2));
ctmatchlatidx(nullidx,1)=0;

% non-canceled trials
nccssdval=sort(unique(nccssd));
nctallmatchlatidx=zeros(length(sacdelay),length(nccssdval));
for ssdval=1:length(nccssdval)
nctallmatchlatidx(:,ssdval)=sacdelay>nccssdval(ssdval)+50 & sacdelay<nccssdval(ssdval)+round(mssrt);
end
% getting ssds for each NNS trial, taking the lowest ssd. 
nctmatchlatidx=zeros(size(nctallmatchlatidx,1),1);
for midx=1:size(nctallmatchlatidx,1)
    if ~isempty(find(nctallmatchlatidx(midx,:),1))
        nctmatchlatidx(midx)=nccssdval(find(nctallmatchlatidx(midx,:),1));
    end
end
end

%% use GUI-independent prealign
getaligndata = prealign(loadfile{:}(1:end-4), trialdirs, tasktype, firstalign,...
     secondalign,  includebad, spikechannel, keepdir,...
     togrey, singlerastplot, [ctmatchlatidx nctmatchlatidx]); % align data, don't plot rasters


%% plots
[allsdf,allrast,allalignidx,allssd,allviscuetimes,allcomp]=disp_cmd([loadfile{:}(1:end-4),'_Clus',num2str(spikechannel)],getaligndata,aligntype,plottype);

%     if strcmp(aligntype,'failed_fast') 
%             %     [p_cancellation,h_cancellation] = cmd_wilco_cancellation(rdd_filename,datalign);
%             disp_cmd([loadfile{:}(1:end-4),'_Clus',num2str(spikechannel)],getaligndata,aligntype,triplot);%0, 0: latmatch, no; triplot, no
%             %     disp_cmd(rdd_filename,datalign,1);
%     elseif strcmp(aligntype,'correct_slow') 
%             disp_cmd([loadfile{:}(1:end-4),'_Clus',num2str(spikechannel)],getaligndata,aligntype,triplot); % keep triplot off until fixed
%     elseif strcmp(aligntype,'ssd') % may need task-specific analysis
% 
%     else
%         
%      end
    catch
        [allsdf,allrast,allalignidx,allssd,allviscuetimes,allcomp]=deal([]);
        continue
    end
    CmdData(FileNb).file=loadfile{:};
    CmdData(FileNb).ssrt=mssrt;
    CmdData(FileNb).tach=[tachomc tachowidth];
    CmdData(FileNb).saclat=sacdelay;
    
    %z-transform the sdf before further analysis by subtracting its nanmean (calculated 
    % from the rasters) and dividing by its standard deviation.
    if strcmp(aligntype,'correct_slow') && size(allsdf,1)>2 
        allsdf=allsdf(1:2,:);
    end
    
    for cp=1:size(allsdf,1)
        for pl=1:size(allsdf,2)
            sessspikes=reshape(allrast{cp,pl}',1,numel(allrast{cp,pl}));
            allsdf{cp,pl}=(allsdf{cp,pl}-nanmean(allsdf{cp,pl}))/nanstd(sessspikes);
        end
    end
    CmdData(FileNb).sdf=allsdf;
    CmdData(FileNb).rast=allrast;
    CmdData(FileNb).algidx=allalignidx;
    CmdData(FileNb).grey=allviscuetimes;
    %get ssd values
    CmdData(FileNb).ssd=allssd;
    
end

%% plotting population
    figure(1); 
    cc=lines(size(CmdData,2)); % one color per file
    if size(cc,1)==8
        cc(8,:)=[0 0.75 0];
    end
   
%first condition plots (for ssd align: canceled vs nss)
    compsdf{1}=nan(size(CmdData,2),size(CmdData(1).sdf{1,1},2));
    compsdf{2}=nan(size(CmdData,2),size(CmdData(1).sdf{2,1},2));
    compsdf{3}=nan(size(CmdData,2),size(CmdData(1).sdf{1,2},2));
    compsdf{4}=nan(size(CmdData,2),size(CmdData(1).sdf{2,2},2));
    
for cmdplots=1:size(CmdData,2)

    for sdfp=1:size(CmdData(cmdplots).sdf,1)
        subplot(size(CmdData(cmdplots).sdf,1),2,sdfp*2-1); hold on 
        plot(CmdData(cmdplots).sdf{sdfp,1},'color',cc(cmdplots,:));
        compsdf{sdfp}(cmdplots,:)=CmdData(cmdplots).sdf{sdfp,1};
    end
end
%legend({CmdData.file},'location','WestOutside');
subplot(size(CmdData(cmdplots).sdf,1),2,2:2:sdfp*2)
compsdf_ci=nanstd(compsdf{1})/ sqrt(size(compsdf{1},1)) ;
patch([1:length(nanmean(compsdf{1})),fliplr(1:length(nanmean(compsdf{1})))],...
    [nanmean(compsdf{1})-compsdf_ci,fliplr(nanmean(compsdf{1})+compsdf_ci)],'b','EdgeColor','none','FaceAlpha',0.1);
hold on;
compsdf_ci=nanstd(compsdf{2})/ sqrt(size(compsdf{2},1)) ;
patch([1:length(nanmean(compsdf{2})),fliplr(1:length(nanmean(compsdf{2})))],...
    [nanmean(compsdf{2})-compsdf_ci,fliplr(nanmean(compsdf{2})+compsdf_ci)],'r','EdgeColor','none','FaceAlpha',0.1);
lineh(1)=plot(nanmean(compsdf{1})); 
lineh(2)=plot(nanmean(compsdf{2}),'r');
currylim=get(gca,'ylim');
patch([plotstart:plotstart+2 fliplr(plotstart:plotstart+2)], ...
            reshape(repmat(currylim,3,1),1,numel(currylim)*3), ...
            [1 0 0],'EdgeColor','none','FaceAlpha',0.5);
%plot SSRT with tach width
for snum=1:size(CmdData,2)
patch([plotstart+CmdData(snum).ssrt:plotstart+CmdData(snum).ssrt+2 fliplr(plotstart+CmdData(snum).ssrt:CmdData(snum).ssrt+plotstart+2 )], ...
            reshape(repmat(currylim,3,1),1,numel(currylim)*3), ...
            [0 0 0],'EdgeColor','none','FaceAlpha',0.5);
end
title('NSS vs CT');
legend(lineh(1:2),'NSS','CT');

%second condition plots (for ssd align: non-canceled vs nss)
 figure(2)
    
for cmdplots=1:size(CmdData,2)

    for sdfp=1:size(CmdData(cmdplots).sdf,1)
        subplot(size(CmdData(cmdplots).sdf,1),2,sdfp*2-1); hold on 
        plot(CmdData(cmdplots).sdf{sdfp,2},'color',cc(cmdplots,:));
        if ~isempty(CmdData(cmdplots).sdf{sdfp,2})
            compsdf{sdfp+size(CmdData(cmdplots).sdf,1)}(cmdplots,:)=CmdData(cmdplots).sdf{sdfp,2};
        else
            continue
        end
   end
end
subplot(size(CmdData(cmdplots).sdf,1),2,2:2:sdfp*2)
compsdf_ci=nanstd(compsdf{1+size(CmdData(cmdplots).sdf,1)})/ sqrt(size(compsdf{1+size(CmdData(cmdplots).sdf,1)},1)) ;
patch([1:length(nanmean(compsdf{1+size(CmdData(cmdplots).sdf,1)})),fliplr(1:length(nanmean(compsdf{1+size(CmdData(cmdplots).sdf,1)})))],...
    [nanmean(compsdf{1+size(CmdData(cmdplots).sdf,1)})-compsdf_ci,fliplr(nanmean(compsdf{1+size(CmdData(cmdplots).sdf,1)})+compsdf_ci)],'b','EdgeColor','none','FaceAlpha',0.1);
hold on;
compsdf_ci=nanstd(compsdf{2+size(CmdData(cmdplots).sdf,1)})/ sqrt(size(compsdf{2+size(CmdData(cmdplots).sdf,1)},1)) ;
patch([1:length(nanmean(compsdf{2+size(CmdData(cmdplots).sdf,1)})),fliplr(1:length(nanmean(compsdf{2+size(CmdData(cmdplots).sdf,1)})))],...
    [nanmean(compsdf{2+size(CmdData(cmdplots).sdf,1)})-compsdf_ci,fliplr(nanmean(compsdf{2+size(CmdData(cmdplots).sdf,1)})+compsdf_ci)],'r','EdgeColor','none','FaceAlpha',0.1);
%put alignment bar
lineh(1)=plot(nanmean(compsdf{1+size(CmdData(cmdplots).sdf,1)})); 
lineh(2)=plot(nanmean(compsdf{2+size(CmdData(cmdplots).sdf,1)}),'r');
currylim=get(gca,'ylim');
patch([plotstart:plotstart+2 fliplr(plotstart:plotstart+2)], ...
            reshape(repmat(currylim,3,1),1,numel(currylim)*3), ...
            [1 0 0],'EdgeColor','none','FaceAlpha',0.5);
%plot SSRT with tach width
for snum=1:size(CmdData,2)
patch([plotstart+CmdData(snum).ssrt:plotstart+CmdData(snum).ssrt+2 fliplr(plotstart+CmdData(snum).ssrt:CmdData(snum).ssrt+plotstart+2 )], ...
            reshape(repmat(currylim,3,1),1,numel(currylim)*3), ...
            [0 0 0],'EdgeColor','none','FaceAlpha',0.5);
end
title('NSS vs NCT');
legend(lineh(1:2),'NSS','NCT');

% diff plots
figure(3); hold on

meandiff(1,:)=nanmean(compsdf{2})-nanmean(compsdf{1});
meandiff(2,:)=nanmean(compsdf{4})-nanmean(compsdf{3});
% first sem
compsdf_ci=nanstd(compsdf{2})/ sqrt(size(compsdf{2},1));     
patch([1:length(nanmean(compsdf{2})),fliplr(1:length(nanmean(compsdf{2})))],...
    [meandiff(1,:)-compsdf_ci,fliplr(meandiff(1,:)+compsdf_ci)],'b','EdgeColor','none','FaceAlpha',0.1);
% second sem
compsdf_ci=nanstd(compsdf{4})/ sqrt(size(compsdf{4},1));
patch([1:length(nanmean(compsdf{4})),fliplr(1:length(nanmean(compsdf{4})))],...
    [meandiff(2,:)-compsdf_ci,fliplr(meandiff(2,:)+compsdf_ci)],'r','EdgeColor','none','FaceAlpha',0.1);
lineh(1)=plot(meandiff(1,:));
lineh(2)=plot(meandiff(2,:),'r');

% looking for statistically significant difference
        diff_pop=meandiff(2,:)-meandiff(1,:);
        diff_preal_epoch=diff_pop(plotstart-200:plotstart);
        difftime_preal=find(abs(diff_preal_epoch)>3*(nanstd(diff_preal_epoch)),1); %three-sigma rule
        diff_postal_epoch=diff_pop(plotstart:plotstart+200);
        difftime_postal=find(abs(diff_postal_epoch)>3*(nanstd(diff_postal_epoch)),1);
            if ~isempty(difftime_preal)
                end_difftime_preal=find(abs(diff_pop(plotstart-difftime_preal+1:end))<=2*(nanstd(diff_preal_epoch)),1);
                end_difftime_preal=plotstart-difftime_preal+end_difftime_preal;
                plot(end_difftime_preal,0,'*b');
                %recursive time search
                difftime_preal=difftime_preal-find(abs(diff_pop(plotstart-200+difftime_preal+1:-1:1))<=2*(nanstd(diff_preal_epoch)),1);
                difftime_preal=plotstart-200+1+difftime_preal;
                plot(difftime_preal,0,'*r');
            end
            if ~isempty(difftime_postal)
                end_difftime_postal=find(abs(diff_pop(plotstart+difftime_postal+1:end))<=2*(nanstd(diff_preal_epoch)),1);
                end_difftime_postal=plotstart+difftime_postal+end_difftime_postal;
                plot(end_difftime_postal,0,'*b');               
                %recursive time search
                difftime_postal=difftime_postal-find(abs(diff_pop(plotstart+difftime_postal+1:-1:1))<=2*(nanstd(diff_preal_epoch)),1);
                difftime_postal=plotstart+1+difftime_postal;
                plot(difftime_postal,0,'*r');
            end
   
currylim=get(gca,'ylim');
% alignment bar
patch([plotstart:plotstart+2 fliplr(plotstart:plotstart+2)], ...
            reshape(repmat(currylim,3,1),1,numel(currylim)*3), ...
            [1 0 0],'EdgeColor','none','FaceAlpha',0.5);
%plot SSRT with tach width
for snum=1:size(CmdData,2)
patch([plotstart+CmdData(snum).ssrt:plotstart+CmdData(snum).ssrt+2 fliplr(plotstart+CmdData(snum).ssrt:CmdData(snum).ssrt+plotstart+2 )], ...
            reshape(repmat(currylim,3,1),1,numel(currylim)*3), ...
            [0 0 0],'EdgeColor','none','FaceAlpha',0.5);
end
title('NC vs NCT');
legend(lineh(1:2),'CT-lmNSS','NCT-lmNSS'); 

% comp='rampNSSvsNCT_ssdalign';
% exportfigname=[cell2mat(regexp(directory,'\w+:\\\w+\\','match')),...
%     'Analysis\Countermanding\',comp];
% %      print(gcf, '-dpng', '-noui', '-opengl','-r600', exportfigname);
% plot2svg([exportfigname,'.svg'],gcf, 'png');

