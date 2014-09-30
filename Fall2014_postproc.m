% keep only filled cells
% (~cellfun('isempty',tasks(:,1)))

%% keep just gapstop data
gsdlist=cellfun(@(x) strcmp(x,'gapstop'),tasks(:,1));
allgsalignmnt=allalignmnt(gsdlist,1);
allgsmssrt=allmssrt(gsdlist,1);
allgspk=allpk(gsdlist,:);
allgsndata=allndata(gsdlist,:); %3 column for 3 aligntype. Each cell has 3 or 4 for diferrent conditions

aligntypes={'failed_fast';'correct_slow','ssd'};

 for gsd=1:size(allgsalignmnt,1)
     
     %ssd alignment
     gsdata=allgsndata(gsd,3);
     if ~isempty(gsdata{:})
         %NSS-CS match latency
         rasters=gsdata{1,1}(1).rasters;
         alignmtt=gsdata{1,1}(1).alignidx;
         start=alignmtt-800; stop=alignmtt+500;
         [sdf, convrasters, convrastsem]=conv_raster(rasters,10,start,stop);
         
         %plots
%          gcf
%          hold on
%          patch([1:length(sdf),fliplr(1:length(sdf))],[sdf-convrastsem,fliplr(sdf+convrastsem)],'k','EdgeColor','none','FaceAlpha',0.1);
%          %plot sdf
%          plot([-800:500],sdf,'Color','b','LineWidth',1.8);
%          hold off

         %normalize by peak activity
         %          cellfun(@(x) x, allgspk(gsd,:),'UniformOutput', false);
         maxpk=max([allgspk{gsd,1}.sac, allgspk{gsd,2}.vis]);
         normsdf=sdf./maxpk;
         
     end
 
 end

