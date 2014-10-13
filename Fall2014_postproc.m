%% keep just gapstop data
gsdlist=cellfun(@(x) strcmp(x,'gapstop'),alltasks(:,1));
allgsalignmnt=allalignmnt(gsdlist,1);
allgsmssrt=allmssrt(gsdlist,1);
allgspk=allpk(gsdlist,:);
allgsndata=allndata(gsdlist,:); %3 column for 3 aligntype. Each cell has 3 or 4 for diferrent conditions

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
    
    for ssdpop=1:2
    popgssdf(ssdpop,:)=nanmean(compgssdf{3}(:,:,ssdpop));
    popgs_ci(ssdpop,:)=std(compgssdf{3}(:,:,ssdpop))/ sqrt(size(compgssdf{3}(:,:,ssdpop),1)) * 1.96;
    
    %% plots
    subplot(1,2,ssdpop)
    patch([1:length(popgssdf(ssdpop,:)),fliplr(1:length(popgssdf(ssdpop,:)))],...
        [popgssdf(ssdpop,:)-popgs_ci(ssdpop,:),fliplr(popgssdf(ssdpop,:)+popgs_ci(ssdpop,:))],'b','EdgeColor','none','FaceAlpha',0.1);
    hold on;
    lineh(ssdpop)=plot(popgssdf(ssdpop,:));
    currylim=get(gca,'ylim');
    patch([800:802 fliplr(800:802)], ...
        reshape(repmat(currylim,3,1),1,numel(currylim)*3), ...
        [1 0 0],'EdgeColor','none','FaceAlpha',0.5);
    end

    
    


    
    
    
    
    
    
    
    
    