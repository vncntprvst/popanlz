listing=dir('B:\Data\Recordings\processed\aligned');
session='S190L6A1'; %'H125L6A2'
alignment='sac'; % ; 'error2sac';'sacstop_cancelstop_non_cancel'
proc='all'; %REX or Sp2
cluster='c1'; % cluster number, e.g. '_c1' 'all'
sesfilesidx=~cellfun('isempty',cellfun(@(x) strfind(x,session), {listing.name},'UniformOutput',false ));
noSHfiles=cellfun('isempty',cellfun(@(x) strfind(x,'2SH'), {listing.name},'UniformOutput',false ));
aligfilesidx=~cellfun('isempty',cellfun(@(x) strfind(x,alignment), {listing.name},'UniformOutput',false ));
if strcmp(proc,'all')
    procfilesidx=1;
else
    procfilesidx=~cellfun('isempty',cellfun(@(x) strfind(x,proc), {listing.name},'UniformOutput',false ));
end
if strcmp(cluster,'all')
    clusfilesidx=1;
else
    clusfilesidx=~cellfun('isempty',cellfun(@(x) strfind(x,cluster), {listing.name},'UniformOutput',false ));
end

FileList={listing(sesfilesidx & noSHfiles & aligfilesidx & procfilesidx & clusfilesidx).name};

% epoch range
allepochrg={[-512 0];[-256 256];[0 512]}; % keep even numbers
epochlabel={'presac';'perisac';'post-sac'};

% cross-correlation window
corrwind=100; % i.e., ±100 ms window

%bins
binwidth=1;  % ms per bin

%prealloc

for filenm=1:length(FileList)
    filename=FileList{filenm};
    
    for epch=1:size(allepochrg,1)
        epochrg=allepochrg{epch};
        
        clear STA STAsem;
        
        % get overall crosscorrelogram, as well as correlation and coherence for
        % x pre-event and y post-event blocks
        try
            [cohrfreq,cohrmag,LFPSumXCorr,Spikenum,ntrials]=SFC(filename,corrwind,epochrg,binwidth);
        catch
            continue;
        end
        
        % STAs
        STA{1,1}=nansum(LFPSumXCorr{1,1})/Spikenum(1); %- mean(nansum(LFPSumXCorr{1,1} )/Spikenum(1));%remove DC offset
        STA{1,2}=nansum(LFPSumXCorr{1,2})/Spikenum(2); %- mean(nansum(LFPSumXCorr{1,2} )/Spikenum(2));%remove DC offset
        STA{2,1}=nanmean(LFPSumXCorr{2,1})/Spikenum(1);% - mean(nanmean(LFPSumXCorr{2,1})/Spikenum(1) );%remove DC offset
        STA{2,2}=nanmean(LFPSumXCorr{2,2})/Spikenum(2);% - mean(nanmean(LFPSumXCorr{2,2})/Spikenum(2) );%remove DC offset
        STA{3,1}=nansum(LFPSumXCorr{3,1})/Spikenum(1); 
        STA{3,2}=nansum(LFPSumXCorr{3,2})/Spikenum(2);
        % STA confidence interval
        STAsem(1,:)=nanstd(LFPSumXCorr{1,1})/ sqrt(size(LFPSumXCorr{1,1},1)); %standard error of the mean
        STAsem(1,:) = STAsem(1,:) * 1.96;
        STAsem(2,:)=nanstd(LFPSumXCorr{1,2})/ sqrt(size(LFPSumXCorr{1,2},1)); %standard error of the mean
        STAsem(2,:) = STAsem(2,:)* 1.96;
                % non-downsampled STA confidence interval
        fullSTAsem(1,:)=nanstd(LFPSumXCorr{3,1})/ sqrt(size(LFPSumXCorr{3,1},1)); %standard error of the mean
        fullSTAsem(1,:) = fullSTAsem(1,:) * 1.96;
        fullSTAsem(2,:)=nanstd(LFPSumXCorr{3,2})/ sqrt(size(LFPSumXCorr{3,2},1)); %standard error of the mean
        fullSTAsem(2,:) = fullSTAsem(2,:)* 1.96;
        
        numcomp=2;
        zcohrmag=cell(2,numcomp);
        zcohrmag_sem=cell(2,numcomp);
        % Fisher z transform
        for compnum=1:numcomp
            zcohrmag{1,compnum}=.5.*log((1+cohrmag{1,compnum})./(1-cohrmag{1,compnum}));
            zcohrmag_sem{1,compnum}=(1/sqrt(ntrials(compnum)-3))*1.96; %95 confidence interval
            zcohrmag{2,compnum}=.5.*log((1+cohrmag{2,compnum})./(1-cohrmag{2,compnum}));
            zcohrmag_sem{2,compnum}=(1/sqrt(ntrials(compnum)-3))*1.96; %95 confidence interval
        end
        
        %plots
        figure('position',[0,0,500,800],'name',[filename '_' epochlabel{epch}],'numbertitle','off')
        %STA
        subplot(2,1,1) 
        patch([1:length(STA{1,1}),fliplr(1:length(STA{1,1}))],[STA{1,1}-STAsem(1,:),fliplr(STA{1,1}+STAsem(1,:))],'b','EdgeColor','none','FaceAlpha',0.1);
        patch([1:length(STA{1,2}),fliplr(1:length(STA{1,2}))],[STA{1,2}-STAsem(2,:),fliplr(STA{1,2}+STAsem(2,:))],'r','EdgeColor','none','FaceAlpha',0.1);
        hold on
        plot(STA{1,1});
%         plot(STA{2,1},'b-.');
        plot(STA{1,2},'r');
%         plot(STA{2,2},'r-.');
        try
            set(gca,'xlim',[1 corrwind*2/binwidth+1],'xtick',1:corrwind:corrwind*2/binwidth+1,'xticklabel',[-corrwind 0 corrwind]); % 'ylim',[nanmin([(STA{1,1}) STA{1,2}]) nanmax([STA{1,1} STA{1,2}])],
        catch
            %bad STA
        end
        title({['STA ' mat2str(corrwind*binwidth) 'ms wd, ' epochlabel{epch}]},'FontSize',16,'FontName','calibri');
        xlabel({'Time (ms)'},'FontSize',16,'FontName','calibri');
        ylabel({'Magnitude'},'FontSize',16,'FontName','calibri');
        legend('','','first condition','','second condition','','location','SouthEastOutside')
        
        % non-downsampled STA
        subplot(2,1,2) 
        patch([1:length(STA{3,1}),fliplr(1:length(STA{3,1}))],[STA{3,1}-fullSTAsem(1,:),fliplr(STA{3,1}+fullSTAsem(1,:))],'b','EdgeColor','none','FaceAlpha',0.1);
        patch([1:length(STA{3,2}),fliplr(1:length(STA{3,2}))],[STA{3,2}-fullSTAsem(2,:),fliplr(STA{3,2}+fullSTAsem(2,:))],'r','EdgeColor','none','FaceAlpha',0.1);
        hold on
        plot(STA{3,1});

        plot(STA{3,2},'r');

        try
            set(gca,'xlim',[1 50*corrwind*2/binwidth+1],'xtick',1:corrwind*50:50*corrwind*2/binwidth+1,'xticklabel',[-corrwind 0 corrwind]); % 'ylim',[nanmin([(STA{1,1}) STA{1,2}]) nanmax([STA{1,1} STA{1,2}])],
        catch
            %bad STA
        end
        title({['STA ' mat2str(corrwind*binwidth) 'ms wd, ' epochlabel{epch}]},'FontSize',16,'FontName','calibri');
        xlabel({'Time (ms)'},'FontSize',16,'FontName','calibri');
        ylabel({'Magnitude'},'FontSize',16,'FontName','calibri');
        legend('','','first condition','','second condition','','location','SouthEastOutside')
        
        % SFC
%         subplot(4,1,3)
%         plot(cohrfreq{2,1}(cohrfreq{2,1}>=0),zcohrmag{2,1}(cohrfreq{2,1}>=0),'b');
%         hold on
%         plot(cohrfreq{2,2}(cohrfreq{2,2}>=0),zcohrmag{2,2}(cohrfreq{2,2}>=0),'r');
%         title({'Coherence (manual method)'},'FontSize',16,'FontName','calibri');
%         xlabel({'Frequency (Hz)'},'FontSize',16,'FontName','calibri');
%         ylabel({'Magnitude'},'FontSize',16,'FontName','calibri');
%                
%         subplot(4,1,4)
%         plot(cohrfreq{1,1},zcohrmag{1,1});
%         hold on
%         plot(cohrfreq{1,2},zcohrmag{1,2},'r');
%         ylim([0 0.5]);  xlim([3 250])
%         title({'z-transformed coherence'},'FontSize',16,'FontName','calibri');
%         xlabel({'Frequency (Hz)'},'FontSize',16,'FontName','calibri');
%         ylabel({'Magnitude'},'FontSize',16,'FontName','calibri');
 
    end
end
