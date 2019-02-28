% Primary use:
% verify that the "flipped" ST/CM profiles are the same neurons
% see data\CommonCells_STCM.xlsx and E:\Dropbox\Vincent Docs\CbTimingPredict\figures\SingleNeuronExample

%% sort and export from Spike2 s (spikes) and t (triggers)
%% plot waveforms
fileName='H146L5A5_23700'; %H146L5A5_23700 H146L5A5_23701
load([fileName 's.mat']);

curVars=who;
loadedData=curVars{cellfun(@(varNames) contains(varNames,fileName),curVars)};
spikeData=eval(loadedData); eval(['clear '  loadedData]);

units=spikeData.codes(:,1);
spiketimes=spikeData.times;
waveforms=spikeData.values;
unitIds=unique(units);unitIds=unitIds(unitIds>0);
legendText=cell(length(unitIds),1);
figure; hold on
for unitNum=1:length(unitIds)
    averageWaveform=mean(waveforms(units==unitIds(unitNum),:));
    plot(averageWaveform,'linewidth',2);
    legendText{unitNum}=['Unit ' num2str(unitNum) ', ' num2str(sum(units==unitIds(unitNum))) ' spikes'];
end
set(gca,'xtick', 0:10:60, 'xticklabel',round(linspace(0,1/30*64,7),1));
xlabel('Time (ms)')
ylabel('Voltage (mV)')
axis('tight');box off;
legend(legendText,'location','southwest')
title(['Average waveforms for ' fileName],'interpreter', 'latex')

wfplot_fileName=[fileName '_AverageWf'];
figDir='E:\Dropbox\Vincent Docs\CbTimingPredict\figures\SingleNeuronExample\';
exportfigname=[figDir wfplot_fileName];
print(gcf, '-dpng', '-noui', '-opengl','-r600', exportfigname);

%% merge file:
REXData=LoadRex_MergeSpk2(fileName);

%% plot rasters for each unit
%% Ecodes
    % ECodes  Countermanding                       Column
%     624y / 427y     Fixation cue                  4
%     644y / 447y     Eye in window                 5
%     664y / 467y     Fixation cue turned off       6
%     684y / 487y     Target onset                  7
%     704y / 507y     Saccade onset or stop signal  8 / 8 or 9, if ecodeout(8)==1503
%     1030            Reward                        11 / 11 (10 if ecodeout(8)~=1503 ?)
%     17385           Bad or aborted trial
%     16386           error code for early saccade  10
%     
     % ECodes  Self-timed                                     Column
%     1001		ENABLECD                                    1
%     602y		Basecode                                    2
%     622y		Onset fix target                            4
%     642y		Eye in window                               5
%     662y		Flash of cue light                          6
%     682y		Cue turned off                              7
%     702y		Rex detected saccade onset                  8
%     722y		Eye is now in target window                 9
%     742y		Re-display target after correct trial       10
%     1030		REWCD                                       11
%     16386     error code for early saccade                8
%     17385     Bad or aborted trial                        9
% Aligning to saccades
% if self-timed saccades, ecode is 702y
% if countermanding task, ecode is 704y
% make rasters
unitIds=unique(vertcat(REXData.Units));
unitIds=unitIds(unitIds>0);
rasters=zeros(size(REXData,2),900,length(unitIds));
for trialNum=1:size(REXData,2)
    findSaccadeEvent=find(floor([REXData(trialNum).Events.Code]/10)==702,1); %704
    if ~isempty(findSaccadeEvent)
        sacTime=REXData(trialNum).Events(findSaccadeEvent).Time-REXData(trialNum).tStartTime; % tStartTime= 1001 code
        epochWindow=REXData(trialNum).SpikeTimes>sacTime-600 & REXData(trialNum).SpikeTimes<sacTime+299;
        for unitNum=1:length(unitIds)
        spikeIndex=REXData(trialNum).Units==unitIds(unitNum);
        spikeTimes=int32(REXData(trialNum).SpikeTimes(epochWindow & spikeIndex))-sacTime+600;
        rasters(trialNum,spikeTimes,unitNum)=1;
        end
    else
        rasters(trialNum,:,:)=NaN;
    end
end
figure; hold on;
for unitNum=1:length(unitIds)
unitRasters=rasters(:,:,unitNum);
unitRasters=unitRasters(~isnan(mean(unitRasters,2)),:);
% figure
% [indy, indx] = ind2sub(size(unitRasters),find(unitRasters)); %find row and column coordinates of spikes
% plot([indx';indx'],[indy';indy'+1],'color','k',...
%     'LineStyle','-','LineWidth',1.8); % plot rasters
unitSDF=fullgauss_filtconv(sum(unitRasters),30,0)./size(unitRasters,1).*1000;
plot(unitSDF)
end
title({['SDF for ' fileName];'Self timed saccade task'},'interpreter', 'latex')
legend(legendText,'location','southwest')

    
    
