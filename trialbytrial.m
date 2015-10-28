function [behav, neur]=trialbytrial(gsdata,conn)

% global directory slash;
% if isempty(directory)
%     [directory,slash]=SetUserDir;
% end

% get cluster indices
[sorted_unit_ids,sunitid_idx]=sort(cellfun(@(x) x.unit_id, gsdata.alldb));
query = ['SELECT c.profile, c.profile_type, c.cluster_id FROM clusters c WHERE cluster_id IN (' sprintf('%.0f,' ,sorted_unit_ids(1:end-1)') num2str(sorted_unit_ids(end)) ')'];
profiles = fetch(conn,query);
sunitid_revidx(sunitid_idx)=1:length(cellfun(@(x) x.unit_id, gsdata.alldb));
neur=mat2cell(([profiles{sunitid_revidx,2}])',ones(110,1));
% clustypes={profiles{sunitid_revidx,1}};
% foo=[profiles{sunitid_revidx,3}];

%% old code
% [peakcct, peaksdf,tbtdircor,tbtdirmsact,tbtdirmsdur]=crosscorel(filename,dataaligned,'all',0);
% tbtcor=nanmedian(abs(tbtdircor))
% tbtcorstd=nanstd(abs(tbtdircor))
% corrcoef([tbtdirmsact tbtdirmsdur])

%% new code
neur=neur(~cellfun('isempty',gsdata.allsacdelay));
neur(:,2)=cellfun(@(x) x(1).rast, gsdata.allndata(~cellfun('isempty',gsdata.allsacdelay), 1),'UniformOutput',false); %NSS trials rast:
neur(:,3)=cellfun(@(x) x(1).alignt, gsdata.allndata(~cellfun('isempty',gsdata.allsacdelay), 1),'UniformOutput',false); %NSS trials alignt:
neur(:,4)=cellfun(@(x) x(1).evttime, gsdata.allndata(~cellfun('isempty',gsdata.allsacdelay), 1),'UniformOutput',false); %NSS trials event time:
%reminder of events: 'cue' 'eyemvt' 'fix' 'rew' 'fail'
behav=cellfun(@(x) x(1).trialnb, gsdata.allndata(~cellfun('isempty',gsdata.allsacdelay), 1),'UniformOutput',false); %NSS trials:
behav(:,2)=cellfun(@(x) [x(2:end).trialnb], gsdata.allndata(~cellfun('isempty',gsdata.allsacdelay), 1),'UniformOutput',false);%all SS trials
behav(:,3)=cellfun(@(x) x.nsst, gsdata.allsacdelay(~cellfun('isempty',gsdata.allsacdelay),1),'UniformOutput',false);%NSS RT
behav(:,4)=cellfun(@(x) x.ncst, gsdata.allsacdelay(~cellfun('isempty',gsdata.allsacdelay),1),'UniformOutput',false);%NCS RT

neur=neur(~cellfun('isempty',behav(:,2)),:);%removing bad apples
behav=behav(~cellfun('isempty',behav(:,2)),:);%removing bad apples

%% behavioral analysis
% impact of previous trial's success/failure on current trial's
% success/failure and on saccade latency

% Vincentizing not necessarily ideal
% see An evaluation of the Vincentizing method of forming group-level response time distributions.
% just pool RT distribution for respective conditions (trials preceded by
% either NSS or SS trial. Anderson-Darling test on distribution http://www.jaqm.ro/issues/volume-6,issue-3/pdfs/1_engmann_cousineau.pdf

behav(:,5)=cellfun(@(x,y,z) x(ismember(y,z+1)), behav(:,3),behav(:,1),behav(:,2),'UniformOutput',false);%RTs for NSS preceded by SS trial
behav(:,6)=cellfun(@(x,y,z) x(ismember(y,z+1)), behav(:,3),behav(:,1),behav(:,1),'UniformOutput',false);%RTs for NSS preceded by NSS trial

[H, adstat, critvalue] = adtest2([behav{:,5}], [behav{:,6}]); %no needt sort RTs
if H==1
    disp(['Distribution are significantly different, A-D stat ' num2str(adstat) ' > critical value ' num2str(critvalue)])
end
figure('Name','Behav_RT_NSS-SSdiff'); ; subplot(1,2,1); hold on;
histh(1)=histogram([behav{:,6}],'Normalization','probability');% mean([behav{:,6}]) %RTs for NSS preceded by NSS trial
histh(2)=histogram([behav{:,5}],'Normalization','probability');% mean([behav{:,5}]) %RTs for NSS preceded by SS trial
histh(1).FaceColor = [0.8 0.2 0.1];
histh(1).EdgeColor = 'w';
histh(2).FaceColor = [0 0.5 0.8];
histh(2).EdgeColor = 'w';
text(300,0.14,{'Two-sample Anderson-Darling'; 'test of significant difference.'; 'alpha = 0.05'})
legend('RTs for NSS preceded by NSS trial','RTs for NSS preceded by SS trial')
legend('boxoff')
xlabel('Reaction times')
title('Normalized distribution histogram')

[N,edges] = histcounts([behav{:,6}], 'Normalization', 'probability');
subplot(1,2,2); plot([edges(1)-diff(edges(1:2))/2 edges(2:end)-diff(edges)/2], cumsum([0 N]))
[N,edges] = histcounts([behav{:,5}], 'Normalization', 'probability');
hold on; plot([edges(1)-diff(edges(1:2))/2 edges(2:end)-diff(edges)/2], cumsum([0 N]))
axis('tight');box off;
xlabel('Reaction times')
title('Cumulative distribution')
legend('RTs for NSS preceded by NSS trial','RTs for NSS preceded by SS trial','location','southeast')
legend('boxoff')

% printing
% cd('E:\BoxSync\Box Sync\Home Folder vp35\Sync\CbTimingPredict\figures')
% exportfigname=get(gcf,'Name');
% exportfigname=strrep(exportfigname,' ','_');
% savefig([exportfigname '.fig']);
% %print png
% print(gcf, '-dpng', '-noui', '-opengl','-r300', exportfigname);
% %print svg
% print(gcf, '-dsvg', '-noui', exportfigname);
% % plot2svg([exportfigname,'.svg'],gcf, 'png');
% close(gcf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Neuron / behavior correlation
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
neur=neur([neur{:,1}]==101,:);
% get slopes (and if really daring, get curve's start as well)
% figure; hold on;
for cell=1:size(neur,1)
    %     figure('Name',['Cell ' num2str(cell)]); hold on;
    for trial=1:size(neur{cell,2},1)
%         figure('Name',['Trial ' num2str(trial)]); hold on;
        % conv % 200 to 500 ms before saccades
        tr_wd_conv = gauss_filtconv(neur{cell,2}(trial,neur{cell,3}-500:neur{cell,3}-1),50);
        %         plot(conv);           
        wd_shift=500-find(tr_wd_conv==max(tr_wd_conv));% plot(tr_wd_conv(1:end-wd_shift));     
        neur{cell,5}(trial)={glmfit(1:500-wd_shift,tr_wd_conv(1:end-wd_shift))};
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

[ccoefs,sigcc]=cellfun(@(x,y) corrcoef(x(2,:)*1000,y),cellfun(@(x) [x{:}], neur(:,5),'UniformOutput',false),behav(:,3),'UniformOutput',false);
ccoefs=[ccoefs{:}];ccoefs=ccoefs(2,1:2:end);
sigcc=[sigcc{:}];sigcc=sigcc(2,1:2:end);
sigcoefs=find(sigcc<0.05); %significant correlation coefficients between pre-sac slope and RT
for cell=1:size(neur,1)
%     neur{sigcoefs(cell),6}=conv_raster(neur{sigcoefs(cell),2}(:,neur{sigcoefs(cell),3}-1900:neur{sigcoefs(cell),3}),20);
%     figure;plot(neur{sigcoefs(sigramp),6});
    [neur{cell,6},~,neur{cell,7}]=conv_raster(neur{cell,2}(:,neur{cell,3}-1900:neur{cell,3}),20);
    % no need for sem actually
end
figure;hold on;
cmap = colormap(gcf);
sigslope_RT=nanmean(cat(1,neur{sigcoefs,6})); sigslope_RT_sem=std(cat(1,neur{sigcoefs,6}))/ sqrt(size(cat(1,neur{sigcoefs,6}),1));
%plot confidence intervals
patch([1:length(sigslope_RT),fliplr(1:length(sigslope_RT))],[sigslope_RT-sigslope_RT_sem,fliplr(sigslope_RT+sigslope_RT_sem)],cmap(1,:),'EdgeColor','none','FaceAlpha',0.1);
%plot sdf
plot(sigslope_RT,'Color',cmap(1,:),'LineWidth',1.8);

nsigslope_RT=nanmean(cat(1,neur{sigcc>=0.05,6})); nsigslope_RT_sem=std(cat(1,neur{sigcc>=0.05,6}))/ sqrt(size(cat(1,neur{sigcc>=0.05,6}),1));
%plot confidence intervals
patch([1:length(nsigslope_RT),fliplr(1:length(nsigslope_RT))],[nsigslope_RT-nsigslope_RT_sem,fliplr(nsigslope_RT+nsigslope_RT_sem)],cmap(20,:),'EdgeColor','none','FaceAlpha',0.1);
%plot sdf
plot(nsigslope_RT,'Color',cmap(20,:),'LineWidth',1.8);


