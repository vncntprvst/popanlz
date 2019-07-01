function [behav, neur]=trialbytrial(cmdata,conn)

% global directory slash;
% if isempty(directory)
%     [directory,slash]=SetUserDir;
% end
plotfig=0;

% get cluster indices
[sorted_unit_ids,sunitid_idx]=sort([cmdata.alldb.unit_id]); %cellfun(@(x) x.unit_id, cmdata.alldb));
query = ['SELECT c.profile, c.profile_type, c.cluster_id FROM clusters c WHERE cluster_id IN (' sprintf('%.0f,' ,sorted_unit_ids(1:end-1)') num2str(sorted_unit_ids(end)) ')'];
profiles = fetch(conn,query);
sunitid_revidx(sunitid_idx)=1:size(cmdata.alldb,1); % length(cellfun(@(x) x.unit_id, cmdata.alldb));%reverse index
% initialize neural data cell array
neur=profiles(sunitid_revidx,2); % much simpler than mat2cell(([profiles{sunitid_revidx,2}])',ones(110,1))  ^^

%% old code
% [peakcct, peaksdf,tbtdircor,tbtdirmsact,tbtdirmsdur]=crosscorel(filename,dataaligned,'all',0);
% tbtcor=nanmedian(abs(tbtdircor))
% tbtcorstd=nanstd(abs(tbtdircor))
% corrcoef([tbtdirmsact tbtdirmsdur])

%% new code
neur=neur(~cellfun('isempty',cmdata.allsacdelay));
neur(:,2,1)=cellfun(@(x) x(1).rast, {cmdata.allndata(~cellfun('isempty',cmdata.allsacdelay)).saccade},...
    'UniformOutput',false); %NSS trials rast, aligned to sac
neur(:,2,2)=cellfun(@(x) x(1).rast, {cmdata.allndata(~cellfun('isempty',cmdata.allsacdelay)).target},...
    'UniformOutput',false); %NSS trials rast, aligned to tgt
neur(:,3,1)=cellfun(@(x) x(1).alignt, {cmdata.allndata(~cellfun('isempty',cmdata.allsacdelay)).saccade},...
    'UniformOutput',false); %NSS trials alignt
neur(:,4,1)=cellfun(@(x) x(1).evttime, {cmdata.allndata(~cellfun('isempty',cmdata.allsacdelay)).saccade},...
    'UniformOutput',false); %NSS trials event time. Reminder of events: 'cue' 'eyemvt' 'fix' 'rew' 'fail'
neur(:,4,2)=cellfun(@(x) x(1).evttime, {cmdata.allndata(~cellfun('isempty',cmdata.allsacdelay)).target},...
    'UniformOutput',false); %NSS trials event time. Reminder of events: 'cue' 'eyemvt' 'fix' 'rew' 'fail'
behav=cellfun(@(x) x(1).trialnb, {cmdata.allndata(~cellfun('isempty',cmdata.allsacdelay)).saccade},...
    'UniformOutput',false)'; %NSS trials:
behav(:,2)=cellfun(@(x) [x(2:end).trialnb], {cmdata.allndata(~cellfun('isempty',cmdata.allsacdelay))...
    .saccade},'UniformOutput',false)';%all SS trials
behav(:,3)=cellfun(@(x) x.nsst, cmdata.allsacdelay(~cellfun('isempty',cmdata.allsacdelay),1),...
    'UniformOutput',false);%NSS RT
behav(:,4)=cellfun(@(x) x.ncst, cmdata.allsacdelay(~cellfun('isempty',cmdata.allsacdelay),1),...
    'UniformOutput',false);%NCS RT
behav(:,5)=cellfun(@(x) x(1).eyevel, {cmdata.allndata(~cellfun('isempty',cmdata.allsacdelay)).target},...
    'UniformOutput',false);%eye velocity
behav(:,8)=cellfun(@(x) x, cmdata.allssds(~cellfun('isempty',cmdata.allsacdelay),1),...
    'UniformOutput',false);%ssds

cellfun(@(x,y) [size(x(2).trialnb,2),size(y{1},1)],...
    {cmdata.allndata(~cellfun('isempty',cmdata.allsacdelay)), cmdata.allssds(~cellfun('isempty',cmdata.allsacdelay),1),...
    'UniformOutput',false);
cellfun(@(x) isempty(x(2).trialnb),{cmdata.allndata(~cellfun('isempty',cmdata.allsacdelay)));

neur=neur(~cellfun('isempty',behav(:,2)),:,:);%removing bad apples
behav=behav(~cellfun('isempty',behav(:,2)),:);%removing bad apples

%% behavioral analysis
% impact of previous trial's success/failure on current trial's
% success/failure and on saccade latency

% Vincentizing not necessarily ideal
% see An evaluation of the Vincentizing method of forming group-level response time distributions.
% just pool RT distribution for respective conditions (trials preceded by
% either NSS or SS trial. Anderson-Darling test on distribution http://www.jaqm.ro/issues/volume-6,issue-3/pdfs/1_engmann_cousineau.pdf

behav(:,6)=cellfun(@(x,y,z) x(ismember(y,z+1)), behav(:,3),behav(:,1),behav(:,2),'UniformOutput',false);%RTs for NSS preceded by SS trial
behav(:,7)=cellfun(@(x,y,z) x(ismember(y,z+1)), behav(:,3),behav(:,1),behav(:,1),'UniformOutput',false);%RTs for NSS preceded by NSS trial
%not going to work:
% behav(:,8)=cellfun(@(x,y,z) x(ismember(z+1,y)), behav(:,8),behav(:,1),behav(:,2),'UniformOutput',false);%SSDs for SS followed by NSS trial

[H, adstat, critvalue] = adtest2([behav{:,6}], [behav{:,7}]); %no needt sort RTs
if H==1
    disp(['Distribution are significantly different, A-D stat ' num2str(adstat) ' > critical value ' num2str(critvalue)])
end
if plotfig
    fh(1)=figure('Name','Behav_RT_NSS-SSdiff'); subplot(1,2,1); hold on;
    histh(1)=histogram([behav{:,7}],'Normalization','probability');% mean([behav{:,7}]) %RTs for NSS preceded by NSS trial
    histh(2)=histogram([behav{:,6}],'Normalization','probability');% mean([behav{:,6}]) %RTs for NSS preceded by SS trial
    histh(1).FaceColor = [0.8 0.2 0.1];
    histh(1).EdgeColor = 'w';
    histh(2).FaceColor = [0 0.5 0.8];
    histh(2).EdgeColor = 'w';
    text(300,0.14,{'Two-sample Anderson-Darling'; 'test of significant difference.'; 'alpha = 0.05'})
    legend('RTs for NSS preceded by NSS trial','RTs for NSS preceded by SS trial')
    legend('boxoff')
    xlabel('Reaction times')
    title('Normalized distribution histogram')
    
    [N,edges] = histcounts([behav{:,7}], 'Normalization', 'probability');
    subplot(1,2,2); plot([edges(1)-diff(edges(1:2))/2 edges(2:end)-diff(edges)/2], cumsum([0 N]))
    [N,edges] = histcounts([behav{:,6}], 'Normalization', 'probability');
    hold on; plot([edges(1)-diff(edges(1:2))/2 edges(2:end)-diff(edges)/2], cumsum([0 N]))
    axis('tight');box off;
    xlabel('Reaction times')
    title('Cumulative distribution')
    legend('RTs for NSS preceded by NSS trial','RTs for NSS preceded by SS trial','location','southeast')
    legend('boxoff')
    
    % make violin plot with gramm
    clear RTdistrib
conditionIdx=cell(length([behav{:,6:7}]),1);
[conditionIdx{sum(~cellfun('isempty', conditionIdx))+1:length([behav{:,7}])}]=deal('RTs for NSS preceded by NSS trial');
[conditionIdx{sum(~cellfun('isempty', conditionIdx))+1:sum(~cellfun('isempty', conditionIdx))+...
    length([behav{:,6}])}]=deal('RTs for NSS preceded by SS trial');

RTdistrib=gramm('x',conditionIdx,'y',[behav{:,7}, behav{:,6}],'color',conditionIdx);
RTdistrib.stat_violin('normalization','area','dodge',0,'fill','edge');%('fill','transparent');
RTdistrib.stat_boxplot('width',0.15);
RTdistrib.set_names('x','Trial conditions','y','RTs');
RTdistrib.set_title('RTs for NSS');
RTdistrib.set_order_options('x',{'RTs for NSS preceded by NSS trial','RTs for NSS preceded by SS trial'});
RTdistrib.set_color_options('map','brewer_dark');
figure;
RTdistrib.draw();
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Neuron / behavior correlation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%try correlation between "slope" and RT
%now with polynomial derivative, but hopefully with nonhomogeneous PP
%% future spiking rate model
% homogeneous poisson coeff
% y=poissrnd(0.5,1,100);
% y(y>0)=1;
% X=1:100;
% coeff = glmfit(X, y,'poisson')
%
% foo=gsdata.allgsndata{1, 1}(1).rast;
% foo=foo(:,200:600);
% [convsdf,convtrial]=conv_raster(foo,10);
% figure;plot(convsdf)
%
% figure; hold on
% for tr=1:size(convtrial,1)
% %     plot(convtrial(43,:))
%     if mean(convtrial(tr,:))>0
%        bla(tr,1:2)= [tr  mean(convtrial(tr,:))];
%     end
% end
%

% coeff = glmfit(1:81, foo(43,180:260),'poisson')
% yfit = glmval(coeff,1:81,-2);
% plot(1:81,yfit,'o')
%% current method
% keep cluster 1 only for the moment
behav=behav([neur{:,1}]==101,:);
neur=neur([neur{:,1}]==101,:,:);
% get slopes (and if really daring, get curve's start as well)
% figure; hold on;
for rec=1:size(neur,1)
    %     figure('Name',['Cell ' num2str(cell)]); hold on;
    for trial=1:size(neur{rec,2},1)
        %         figure('Name',['Trial ' num2str(trial)]); hold on;
        % conv % 200 to 500 ms before saccades
        tr_wd_conv = gauss_filtconv(neur{rec,2,1}(trial,neur{rec,3,1}-500:neur{rec,3,1}-1),50);
        %         plot(conv);
        wd_shift=500-find(tr_wd_conv==max(tr_wd_conv));% plot(tr_wd_conv(1:end-wd_shift));
        neur{rec,5,1}(trial)={glmfit(1:500-wd_shift,tr_wd_conv(1:end-wd_shift))};
        %         yfit = glmval(coeff,1:301-wd_shift,'identity');
        
        %         pf = polyfit(1:301, neur{cell,2}(trial,neur{cell,3}-500:neur{cell,3}-200),4);
        %         k = polyder(pf);plot(1:301,(1:301)*k(end-1)+k(end));
        %         coeff = glmfit(1:301,neur{cell,2}(trial,neur{cell,3}-500:neur{cell,3}-200),'poisson')
        %         yfit = glmval(coeff,1:301,'log');
        %         plot(yfit);
        
    end
end
% figure; histogram(cellfun(@(x) x(2)*1000,neur(:,5)),10)
%reminder: neur 'labels' are 'cluster id','NSST raster','NSST alignment','NSST evt time','slopes'

[ccoefs,sigcc]=cellfun(@(x,y) corrcoef(x(2,:),y),cellfun(@(x) [x{:}], neur(:,5,1),'UniformOutput',false),behav(:,3),'UniformOutput',false);
ccoefs=[ccoefs{:}];ccoefs=ccoefs(2,1:2:end);
sigcc=[sigcc{:}];sigcc=sigcc(2,1:2:end);
sigcoefs=sigcc<0.05; %significant correlation coefficients between pre-sac slope and RT
disp(['mean significant coeff ' num2str(mean(ccoefs(sigcoefs))) ' +/- ' num2str(std(ccoefs(sigcoefs))/sqrt(size(ccoefs(sigcoefs),2)))])
sigma=20;
for rec=1:size(neur,1)
    %     neur{sigcoefs(cell),6}=conv_raster(neur{sigcoefs(cell),2}(:,neur{sigcoefs(cell),3}-1900:neur{sigcoefs(cell),3}),20);
    %     figure;plot(neur{sigcoefs(sigramp),6});
    [neur{rec,6,1},~,neur{rec,7,1}]=conv_raster(neur{rec,2,1}(:,neur{rec,3,1}-1900:neur{rec,3,1}+3*sigma),sigma);
    % no need for sem actually
end
if plotfig
    fh(2)=figure('Name','ramp_RT correlation - Mean pre-sac activity');hold on;
    cmap = colormap(gcf);
    sigslope_RT=nanmean(cat(1,neur{sigcoefs,6,1})); sigslope_RT_sem=std(cat(1,neur{sigcoefs,6,1}))/ sqrt(size(cat(1,neur{sigcoefs,6,1}),1));
    %plot confidence intervals
    patch([1:length(sigslope_RT),fliplr(1:length(sigslope_RT))],[sigslope_RT-sigslope_RT_sem,fliplr(sigslope_RT+sigslope_RT_sem)],cmap(1,:),'EdgeColor','none','FaceAlpha',0.1);
    %plot sdf
    ploth(1)=plot(sigslope_RT,'Color',cmap(1,:),'LineWidth',2,'DisplayName','significant RT/slope correlation');
    
    nsigslope_RT=nanmean(cat(1,neur{sigcc>=0.05,6,1})); nsigslope_RT_sem=std(cat(1,neur{sigcc>=0.05,6,1}))/ sqrt(size(cat(1,neur{sigcc>=0.05,6,1}),1));
    %plot confidence intervals
    patch([1:length(nsigslope_RT),fliplr(1:length(nsigslope_RT))],[nsigslope_RT-nsigslope_RT_sem,fliplr(nsigslope_RT+nsigslope_RT_sem)],cmap(20,:),'EdgeColor','none','FaceAlpha',0.1);
    %plot sdf
    ploth(2)=plot(nsigslope_RT,'Color',cmap(20,:),'LineWidth',2,'DisplayName','non-significant RT/slope correlation');
    
    [~,ks_p] = kstest2(sigslope_RT,nsigslope_RT); % test should focus on peak times, rather than sdf over whole time course
    text(10,max(get(gca,'ylim'))-10,{'Two-sample Kolmogorov-Smirnov test';['p=' num2str(ks_p)]})
    
    title('population average pre-sac sdf')
    xlabel('Time')
    ylabel('Firing rate')
    set(gca,'xticklabel',-(max(get(gca,'xlim'))):200:0,'TickDir','out'); %'xtick',1:100:max ...
    box off; % axis('tight');
    legh=legend(ploth(1:2),'location','southeast');
end

%% Trial-by-Trial Correlation of Neural and Behavioral Latency (a la Erlich)
[pkOffsets,recIdx]=latency_analysis(behav,neur,0);
% test correlation between FR and eye mvt
[ccoefs,sigcc]=cellfun(@(x,y) corrcoef(x,y),pkOffsets(:,1),pkOffsets(:,2),'UniformOutput',false);
ccoefs=[ccoefs{:}];ccoefs=ccoefs(2,1:2:end);[ccoefs(isnan(ccoefs))]=deal(0);
sigcc=[sigcc{:}];sigcc=sigcc(2,1:2:end);[sigcc(isnan(sigcc))]=deal(1);
sigcoefs=sigcc<0.05; %significant trial-by-trial correlation of neural and behavioral latency
if plotfig
%     figure; hold on;
% histogram(ccoefs,-1:0.1:1);
% histogram(ccoefs(sigcoefs),-1:0.1:1);
% % patch([linspace(mean(ccoefs)-0.1,mean(ccoefs)+0.1,length(min(get(gca,'ylim')):max(get(gca,'ylim')))),...
% %     fliplr(linspace(mean(ccoefs)-0.1,mean(ccoefs)+0.1,length(min(get(gca,'ylim')):max(get(gca,'ylim')))))],...
% %     [min(get(gca,'ylim')):max(get(gca,'ylim')),fliplr(min(get(gca,'ylim')):max(get(gca,'ylim')))],cmap(1,:),'EdgeColor','none','FaceAlpha',1);
% plot(mean(ccoefs),max(get(gca,'ylim'))/1.5,'kd')
    recIdx=find(recIdx);
    fh(3)=figure('Name','FR - Eye velocity correlation');hold on;
    cmap = colormap(gcf);
    %% stoped here. Replace sigslope_RT with FReyevel_cor
    sigslope_RT=nanmean(cat(1,neur{recIdx(sigcoefs),6,1})); sigslope_RT_sem=std(cat(1,neur{recIdx(sigcoefs),6,1}))/ sqrt(size(cat(1,neur{recIdx(sigcoefs),6,1}),1));
    %plot confidence intervals
    patch([1:length(sigslope_RT),fliplr(1:length(sigslope_RT))],[sigslope_RT-sigslope_RT_sem,fliplr(sigslope_RT+sigslope_RT_sem)],cmap(1,:),'EdgeColor','none','FaceAlpha',0.1);
    %plot sdf
    ploth(1)=plot(sigslope_RT,'Color',cmap(1,:),'LineWidth',2,'DisplayName','significant RT/slope correlation');
    
    nsigslope_RT=nanmean(cat(1,neur{sigcc>=0.05,6,1})); nsigslope_RT_sem=std(cat(1,neur{sigcc>=0.05,6,1}))/ sqrt(size(cat(1,neur{sigcc>=0.05,6,1}),1));
    %plot confidence intervals
    patch([1:length(nsigslope_RT),fliplr(1:length(nsigslope_RT))],[nsigslope_RT-nsigslope_RT_sem,fliplr(nsigslope_RT+nsigslope_RT_sem)],cmap(20,:),'EdgeColor','none','FaceAlpha',0.1);
    %plot sdf
    ploth(2)=plot(nsigslope_RT,'Color',cmap(20,:),'LineWidth',2,'DisplayName','non-significant RT/slope correlation');
    
    [~,ks_p] = kstest2(sigslope_RT,nsigslope_RT); % test should focus on peak times, rather than sdf over whole time course
    text(10,max(get(gca,'ylim'))-10,{'Two-sample Kolmogorov-Smirnov test';['p=' num2str(ks_p)]})
    
    title('population average pre-sac sdf')
    xlabel('Time')
    ylabel('Firing rate')
    set(gca,'xticklabel',-(max(get(gca,'xlim'))):200:0,'TickDir','out'); %'xtick',1:100:max ...
    box off; % axis('tight');
    legh=legend(ploth(1:2),'location','southeast');
end

% test correlation between FR and stop signal delay from previous SST
[ccoefs,sigcc]=cellfun(@(x,y) corrcoef(x,y),pkOffsets(recIdx,3),pkOffsets(recIdx,4),'UniformOutput',false);
gCoeffs=cellfun(@(x) size(x,1)>1, ccoefs);
ccoefs=ccoefs(gCoeffs);sigcc=sigcc(gCoeffs);
ccoefs=[ccoefs{:}];ccoefs=ccoefs(2,1:2:end);[ccoefs(isnan(ccoefs))]=deal(0);
sigcc=[sigcc{:}];sigcc=sigcc(2,1:2:end);[sigcc(isnan(sigcc))]=deal(1);
sigcoefs=sigcc<0.05; %significant trial-by-trial correlation of neural and behavioral latency
figure; hold on;
histogram(ccoefs,-1:0.1:1);
histogram(ccoefs(sigcoefs),-1:0.1:1);
% patch([linspace(mean(ccoefs)-0.1,mean(ccoefs)+0.1,length(min(get(gca,'ylim')):max(get(gca,'ylim')))),...
%     fliplr(linspace(mean(ccoefs)-0.1,mean(ccoefs)+0.1,length(min(get(gca,'ylim')):max(get(gca,'ylim')))))],...
%     [min(get(gca,'ylim')):max(get(gca,'ylim')),fliplr(min(get(gca,'ylim')):max(get(gca,'ylim')))],cmap(1,:),'EdgeColor','none','FaceAlpha',1);
plot(mean(ccoefs),max(get(gca,'ylim'))/1.5,'kd')

%test pooling all offsets together
pooledPkOffsets=vertcat(pkOffsets{:,3});
pooledPkOffsets(:,2)=vertcat(pkOffsets{:,4});
[ccoefs,sigcc]=corrcoef(foo);

% printing
cd('E:\BoxSync\Box Sync\Home Folder vp35\Sync\CbTimingPredict\figures')
exportfigname=get(fh(2),'Name');
exportfigname=strrep(exportfigname,' ','_');
savefig([exportfigname '.fig']);
%print png
print(gcf, '-dpng', '-noui', '-opengl','-r300', exportfigname);
%print svg
print(gcf, '-dsvg', '-noui', exportfigname);
% plot2svg([exportfigname,'.svg'],gcf, 'png');
close(gcf)

%% Single trial / trial-by trial examples
% use data from 1 cluster, and maybe non-significant RT/slope correlation
% make trial progression RT vs peak figure, if one can segment session in
% different period (see staircase delays)
% make scatter plot with distribution of peak vs RT, following NSST / SST
% trial-by-trial variability: 
%     Krakauer: https://www.physiology.org/doi/10.1152/jn.00897.2013
% reduction in trial-by-trial variability for the trained task 
% but not for an untrained control task, without changes in mean trajectory. 
% The decrease in variability for the trained task was associated with increased 
% activation in contralateral primary motor and premotor cortical areas and 
% in ipsilateral cerebellum.
%  -> do RT or trajectory variablity stabilize, and so response properties with it?  
% Anne Churchland:  http://www.jneurosci.org/content/37/19/4954
% trial-to-trial spike count variability during stimulus presentation and decision formation
% Variability decreased sharply, suggesting that the network is stabilized by inputs
% -> are dentate neurons activity during sensory and delay processing stabilized by inputs?
% Chris Gomez:  evidence for the trial-by-trial strategy was obtained using 
% random-effects GLM models with the speed or angle change category and group
% as independent factors and the deviation from the optimal point as the dependent variable

% Decoding
% See https://www.nature.com/articles/nrn2578
% Single-neuron discrimination analysis: https://www.nature.com/articles/nature12236

% Could use ROC analysis like Robin (or Chen 2015) to find out time of divergence. But
% this is not a discrimination task. Maybe do as follow:
%  Chris Harvey: https://www.nature.com/articles/nature10918
% To quantify the trajectory divergence, we used a classifier
% based on the distance from an individual trial trajectory to the mean right 
% and left choice trajectories at single time points
% -> Neuronal circuit trajectories on correct and erroneous SST 
% Use classifier to find critical period at which decision is made


