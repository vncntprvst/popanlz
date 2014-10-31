%% keep just gapstop data
gsdlist=cellfun(@(x) strcmp(x,'gapstop'),alltasks(:,1));
allgsalignmnt=allalignmnt(gsdlist,1);
allgsmssrt=allmssrt(gsdlist,1);
allgspk=allpk(gsdlist,:);
allgsndata=allndata(gsdlist,:); %3 column for 3 aligntype. Each cell has 3 or 4 for diferrent conditions
allgs_rec_id=all_rec_id(gsdlist,1);
allgsstats=allstats(gsdlist,1);

%% cluster analysis of population
%convolve rasters with 200ms before saccade, 200 after saccade, 10ms kernel
%time window (so add 30 ms at both ends, which will be cut)

sacresps=cellfun(@(x) conv_raster(x(1,1).rast,10,x(1,1).alignt-230,x(1,1).alignt+229), allgsndata(:,1), 'UniformOutput',false);
sacresps=cat(1,sacresps{:});

seeds=cellfun(@(x) mean(x(1,100:200))-mean(x(1,200:300)), mat2cell(sacresps,ones(size(sacresps,1),1)));
wavedropseed=seeds==max(seeds);
waveburstseed=seeds==min(seeds);
waveflatseed=seeds==min(abs(seeds));
seeds=[sacresps(wavedropseed,:);...
     sacresps(waveburstseed,:)];
%     sacresps(waveflatseed,:)];

% %plot seeds
figure
plot(seeds(1,:))
hold on
plot(seeds(2,:),'r')
% plot(seeds(3,:),'g')
    
% at random
% [IDX,C,sumd,D]=kmeans(sacresps,3,'dist','city','display','iter'); 
% seeded
%[IDX,C,sumd,D]=kmeans(sacresps,2,'dist','city','start',seeds,'display','iter'); 
%plot means
figure
plot(C(1,:),'b')
hold on
plot(C(2,:),'r')
% plot(C(3,:),'g')

%plot means
figure
plot(C(1,:),'b')
hold on
plot(C(2,:),'r')
% plot(C(3,:),'g')
%scatterplot
figure
scatter3(D(IDX==1,1),D(IDX==1,2),D(IDX==1,3),40,'b.')
hold on
scatter3(D(IDX==2,1),D(IDX==2,2),D(IDX==2,3),40,'r.')
scatter3(D(IDX==3,1),D(IDX==3,2),D(IDX==3,3),40,'g.')

% %plot sdfs from each cluster
clus1=find(IDX==1); clus2=find(IDX==2); clus3=find(IDX==3);
figure('name','cluster1')
subplotdim=[ceil(length(clus1)/2)-(2*floor(length(clus1)/10)),2+floor(length(clus1)/10)];
for sacplot=1:length(clus1)
    subplot(subplotdim(1),subplotdim(2),sacplot)
    plot(sacresps(clus1(sacplot),:));
    ylim=get(gca,'ylim');
    set(gca,'ylim',[0 ylim(2)]);
    text(20,10,['sacplot ' num2str(clus1(sacplot))]);
end
figure('name','cluster2')
subplotdim=[ceil(length(clus2)/2)-(2*floor(length(clus2)/10)),2+floor(length(clus2)/10)];
for sacplot=1:length(clus2)
    subplot(subplotdim(1),subplotdim(2),sacplot)
    plot(sacresps(clus2(sacplot),:));
    ylim=get(gca,'ylim');
    set(gca,'ylim',[0 ylim(2)]);
    text(20,10,['sacplot ' num2str(clus2(sacplot))]);
end
figure('name','cluster3')
subplotdim=[ceil(length(clus3)/2)-(2*floor(length(clus3)/10)),2+floor(length(clus3)/10)];
for sacplot=1:length(clus3)
    subplot(subplotdim(1),subplotdim(2),sacplot)
    plot(sacresps(clus3(sacplot),:));
    ylim=get(gca,'ylim');
    set(gca,'ylim',[0 ylim(2)]);
    text(20,10,['sacplot ' num2str(clus3(sacplot))]);
end

%PCA
[coeffs,score,latent] = pca(sacresps);
% D=coeffs(:,1:8)';
figure
plot(coeffs(:,1),'.','MarkerSize',.5)
hold on 
plot(coeffs(:,2),'.','MarkerSize',.5)
plot(coeffs(:,3),'.','MarkerSize',.5)
plot(coeffs(:,4),'.','MarkerSize',.5)

%scatterplot
% in 2D
figure
plot(score(:,1),score(:,2),'.b')
hold on
scatter3(score([6,11,21,24,29,32],1),score([6,11,21,24,29,32],2),score([6,11,21,24,29,32],3),80,'r','filled')
scatter3(score([7,14,17,25,27],1),score([7,14,17,25,27],2),score([7,14,17,25,27],3),80,'g','filled')

% in 2D, with variance multiplier
figure
plot(score(:,1).*latent(1),score(:,2).*latent(2),'.b')

% in 3D
figure
scatter3(score(:,1),score(:,2),score(:,3),40,'b.')
hold on
% outline just seeds
scatter3(score(wavedropseed,2),score(wavedropseed,3),score(wavedropseed,4),80,'r','filled')
scatter3(score(waveburstseed,2),score(waveburstseed,3),score(waveburstseed,4),80,'g','filled')
scatter3(score(waveflatseed,2),score(waveflatseed,3),score(waveflatseed,4),80,'k','filled')
% outline bests
scatter3(score([6,11,21,24,29,32],1),score([6,11,21,24,29,32],2),score([6,11,21,24,29,32],3),80,'r','filled')
scatter3(score([7,14,17,25,27],1),score([7,14,17,25,27],2),score([7,14,17,25,27],3),80,'g','filled')


%
scatter3(score(IDX==1,1),score(IDX==1,2),score(IDX==1,3),40,'b.')
hold on
scatter3(score(IDX==2,1),score(IDX==2,2),score(IDX==2,3),40,'r.')
scatter3(score(IDX==3,1),score(IDX==3,2),score(IDX==3,3),40,'g.')

    %% colors for population plots
%     figure(1);
    cc=lines(size(allgsalignmnt,1)); % one color per file
    if size(cc,1)==8
        cc(8,:)=[0 0.75 0];
    end
    
    %% prealloc compile data
    compgssdf{1}=nan(size(allgsalignmnt,1),1301,3);
    compgssdf{2}=nan(size(allgsalignmnt,1),901,3);
    compgssdf{3}=nan(size(allgsalignmnt,1),1301,4);
    compgssdf{4}=nan(size(allgsalignmnt,1),1001,3);
    compgssdf{5}=nan(size(allgsalignmnt,1),1001,3);
    
    %% get sdf's
 for gsd=1:size(allgsalignmnt,1) %compile data across all files for each condition

     %% sac alignment ('failed_fast')
     gsdata=allgsndata{gsd,1}; %this will be 3 structures containing rasters 
                               % and alignments sac/cancellation/ wrong sac)
     try 
         gspk=allgspk{gsd,1}.sac;
     catch nopeak
         gspk=0;
     end
     
     if ~isempty(gsdata) && gspk~=0
         for sacalg=1:3
             try
             rasters=gsdata(sacalg).rast;
             alignmtt=gsdata(sacalg).alignt;
             start=alignmtt-900; stop=alignmtt+400;
             [sdf, convrasters, convrastsem]=conv_raster(rasters,10,start,stop); 
             
             
             %normalize sdf by peak activity
             normsdf=sdf./gspk;

             %% plots
    %          figure(1)
    %          hold off
    %          patch([1:length(sdf),fliplr(1:length(sdf))],[sdf-convrastsem,fliplr(sdf+convrastsem)],'k','EdgeColor','none','FaceAlpha',0.1);
    %          hold on
    %          %plot sdf
    %          plot(sdf,'Color','b','LineWidth',1.8);
    %          set(gca,'xtick',[1:100:1301],'xticklabel',[-800:100:500])
    %          close(gcf)

             %% store 
             compgssdf{1}(gsd,:,sacalg)=sdf; 
             catch norast
             end
         end
     else
         %keep nans
     end
     
     %% tgt alignment ('correct_slow')
     gsdata=allgsndata{gsd,2}; %this will be 3 structures containing rasters 
                               % and alignments tgt/tgt-CS/tgt-NCS)
     try 
         gspk=allgspk{gsd,2}.vis;
     catch nopeak
         gspk=0;
     end
     
     if ~isempty(gsdata) && gspk~=0
         for tgtalg=1:3
             try
             rasters=gsdata(tgtalg).rast;
             alignmtt=gsdata(tgtalg).alignt;
             start=alignmtt-200; stop=alignmtt+700;
             [sdf, convrasters, convrastsem]=conv_raster(rasters,10,start,stop);

             %normalize sdf by peak activity
             normsdf=sdf./gspk;

             %% plots
    %          figure(1)
    %          hold off
    %          patch([1:length(sdf),fliplr(1:length(sdf))],[sdf-convrastsem,fliplr(sdf+convrastsem)],'k','EdgeColor','none','FaceAlpha',0.1);
    %          hold on
    %          %plot sdf
    %          plot(sdf,'Color','b','LineWidth',1.8);
    %          set(gca,'xtick',[1:100:1301],'xticklabel',[-800:100:500])
    %          close(gcf)

             %% store 
             compgssdf{2}(gsd,:,tgtalg)=sdf; 
             catch norast
             end
         end
     else
         %keep nans
     end
     
     %% ssd alignment ('ssd')
     gsdata=allgsndata{gsd,3}; %this will be 4 structures containing rasters 
                               % and alignment CS & NCS with corresponding 
                               % Lm-NSS trials (Lm-CS,Lm-NCS,CS,NCS)
     try 
         gspk=allgspk{gsd,3}.vis;
     catch nopeak
         gspk=0;
     end
     
     if ~isempty(gsdata) && gspk~=0
         for ssdalg=1:4
             try
             rasters=gsdata(ssdalg).rast;
             alignmtt=gsdata(ssdalg).alignt;
             start=alignmtt-800; stop=alignmtt+500;
             [sdf, convrasters, convrastsem]=conv_raster(rasters,10,start,stop);

             %normalize sdf by peak activity
             normsdf=sdf./gspk;

             %% plots
    %          figure(1)
    %          hold off
    %          patch([1:length(sdf),fliplr(1:length(sdf))],[sdf-convrastsem,fliplr(sdf+convrastsem)],'k','EdgeColor','none','FaceAlpha',0.1);
    %          hold on
    %          %plot sdf
    %          plot(sdf,'Color','b','LineWidth',1.8);
    %          set(gca,'xtick',[1:100:1301],'xticklabel',[-800:100:500])
    %          close(gcf)

             %% store 
             compgssdf{3}(gsd,:,ssdalg)=sdf; 
             catch norast
             end
         end
     else
         %keep nans
     end
 
     %% corrective saccade alignment ('corrsacfailed')
     gsdata=allgsndata{gsd,4}; %this will be 3 structures containing rasters 
                               % and alignments corsac/~corsac/NCS corsac)
     try 
         gspk=allgspk{gsd,4}.corsac;
     catch nopeak
         gspk=0;
     end
     
     if ~isempty(gsdata) && gspk~=0
         for csacalg=1:3
             try
             rasters=gsdata(csacalg).rast;
             alignmtt=gsdata(csacalg).alignt;
             start=alignmtt-800; stop=alignmtt+200;
             [sdf, convrasters, convrastsem]=conv_raster(rasters,10,start,stop);

             %normalize sdf by peak activity
             normsdf=sdf./gspk;

             %% plots
    %          figure(1)
    %          hold off
    %          patch([1:length(sdf),fliplr(1:length(sdf))],[sdf-convrastsem,fliplr(sdf+convrastsem)],'k','EdgeColor','none','FaceAlpha',0.1);
    %          hold on
    %          %plot sdf
    %          plot(sdf,'Color','b','LineWidth',1.8);
    %          set(gca,'xtick',[1:100:1301],'xticklabel',[-800:100:500])
    %          close(gcf)

             %% store 
             compgssdf{4}(gsd,:,csacalg)=sdf; 
             catch norast
             end
         end
     else
         %keep nans
     end
     
     %% reward alignment ('rewcorrect_rewslow')
     gsdata=allgsndata{gsd,5}; %this will be 3 structures containing rasters 
                               % and alignments corsac/~corsac/NCS corsac)
     try 
         gspk=allgspk{gsd,5}.rew;
     catch nopeak
         gspk=0;
     end
     
     if ~isempty(gsdata) && gspk~=0
         for rewalg=1:3
             try
             rasters=gsdata(rewalg).rast;
             alignmtt=gsdata(rewalg).alignt;
             start=alignmtt-800; stop=alignmtt+200;
             [sdf, convrasters, convrastsem]=conv_raster(rasters,10,start,stop);

             %normalize sdf by peak activity
             normsdf=sdf./gspk;

             %% plots
    %          figure(1)
    %          hold off
    %          patch([1:length(sdf),fliplr(1:length(sdf))],[sdf-convrastsem,fliplr(sdf+convrastsem)],'k','EdgeColor','none','FaceAlpha',0.1);
    %          hold on
    %          %plot sdf
    %          plot(sdf,'Color','b','LineWidth',1.8);
    %          set(gca,'xtick',[1:100:1301],'xticklabel',[-800:100:500])
    %          close(gcf)

             %% store 
             compgssdf{5}(gsd,:,rewalg)=sdf; 
             catch norast
             end
         end
     else
         %keep nans
     end
 end

    %% population activity, ci and plots
    
    ssdfig=figure('name','ssd figure');
    
    popgssdf=NaN(4,size(compgssdf{3},2));
    popgs_ci=NaN(4,size(compgssdf{3},2));
    
    %% ssd pop
    %1st plots
    subplot(1,2,1)
    
    for ssdpop=1:2:3
        popgssdf(ssdpop,:)=nanmean(compgssdf{3}(:,:,ssdpop));
    popgs_ci(ssdpop,:)=nanstd(compgssdf{3}(:,:,ssdpop))/ sqrt(size(compgssdf{3}(:,:,ssdpop),1));

    patch([1:length(popgssdf(ssdpop,:)),fliplr(1:length(popgssdf(ssdpop,:)))],...
        [popgssdf(ssdpop,:)-popgs_ci(ssdpop,:),fliplr(popgssdf(ssdpop,:)+popgs_ci(ssdpop,:))],...
        cc(ssdpop,:),'EdgeColor','none','FaceAlpha',0.1);
    hold on;
    lineh(ssdpop)=plot(popgssdf(ssdpop,:),'color',cc(ssdpop,:));
    end
    currylim=get(gca,'ylim');
    patch([798:802 fliplr(798:802)], ...
        reshape(repmat([0 currylim(2)],5,1),1,numel(currylim)*5), ...
        [0 0 0],'EdgeColor','none','FaceAlpha',0.5);
    
    %2nd plots
    subplot(1,2,2)
    
    for ssdpop=2:2:4
    popgssdf(ssdpop,:)=nanmean(compgssdf{3}(:,:,ssdpop));
    popgs_ci(ssdpop,:)=nanstd(compgssdf{3}(:,:,ssdpop))/ sqrt(size(compgssdf{3}(:,:,ssdpop),1));

    patch([1:length(popgssdf(ssdpop,:)),fliplr(1:length(popgssdf(ssdpop,:)))],...
        [popgssdf(ssdpop,:)-popgs_ci(ssdpop,:),fliplr(popgssdf(ssdpop,:)+popgs_ci(ssdpop,:))],...
        cc(ssdpop,:),'EdgeColor','none','FaceAlpha',0.1);
    hold on;
    lineh(ssdpop)=plot(popgssdf(ssdpop,:),'color',cc(ssdpop,:));
    end
        currylim=get(gca,'ylim');
    patch([798:802 fliplr(798:802)], ...
        reshape(repmat([0 currylim(2)],5,1),1,numel(currylim)*5), ...
        [0 0 0],'EdgeColor','none','FaceAlpha',0.7);

    
    


    
    
    
    
    
    
    
    
    