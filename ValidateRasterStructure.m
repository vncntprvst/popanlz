%% definitions
clearvars; clear global;
taskType= 'st_saccades'; %gapstop
% ecodes for coumtermanding task
[fixcode, fixoffcode, tgtcode, tgtoffcode, saccode, ...
    stopcode, rewcode, ~, errcode1, errcode2, errcode3, basecode, ...
    dectgtcode, decsaccode] = taskfindecode(taskType); 
% processed files directory
procDir='E:\Data\Recordings\processed\';
% alignment types
alignTypes={'saccade';'target';'reward'}; %'corrective saccade' 'stopsignal'
allProcData=struct('saccade',struct([]),'target',struct([]),...
    'reward',struct([])); % 'stopsignal',struct([]),
%% load cmdata
load('E:\Dropbox\Vincent Docs\CbTimingPredict\data\cDn_stdata.mat'); % top_cortex_cmdata %cDn_cmdata
if ~isfield(stdata,'fileNames')
    try
        stdata.fileNames=stdata.alldb.recNames;
    catch
        listRecIds=cellfun(@(x) x.rec_id,stdata.alldb);
        CCNdb = connect2DB('vp_sldata');
        query = ['SELECT r.a_file FROM recordings r WHERE recording_id IN (' ...
    sprintf('%.0f,' ,cellfun(@(x) x.rec_id,stdata.alldb(1:end-1,1))) num2str(stdata.alldb{end,1}.rec_id) ')'];
%         query = ['select r.e_file from recordings r where recording_id=' ...
%             num2str(stData.alldb{commonCells.st(commonFileNum), 1}.rec_id) ';'];
        st_fileNames=fetch(CCNdb,query);
        [~,resort]=sort(listRecIds);[~,resort]=sort(resort);st_fileNames=st_fileNames(resort); %already in ascending order, but just to make sure
        stdata.fileNames=st_fileNames; %table(st_fileNames,'VariableNames',{'fileNames'});
    end
end
    %% go file by file
    numFiles=size(stdata.fileNames  ,1);
    for fileNum=1:numFiles
        % get data for that file from the structure
        structData = stdata.allndata(fileNum,:); %,[1,2,3,5]); %cmdata.allndata(fileNum);
        fileName=stdata.fileNames{fileNum}(1:end-1);
        switch fileName(1)
            case 'R'
                dirName='Rigel';
            case 'S'
                dirName='Sixx';
            case 'H'
                dirName='Hilda';
        end
        %load corresponding processed file
        try
            fileDir=[procDir dirName];
            load([fileDir filesep fileName '_REX.mat']);
        catch
            %select file to load
            [fileName,fileDir] = uigetfile({'*.mat','Processed Files';...
                '*.*','All Files' },'Select file to load',...
                [fileDir filesep fileName '_REX.mat']);
        end
        %% now go alignment by alignment 
        for alignTypeNum=1:3 % saccade / target / reward / %1:4 (columns : saccade / target / stop signal / reward / ) % corrective saccade
            % go trial type by trial type
            try 
                numTrialTypes= size(structData{1,alignTypeNum},2); %size(structData.(alignTypes{alignTypeNum}),2);
            catch 
                numTrialTypes=1;
            end
            if strcmp(taskType,'gapstop')
                if numTrialTypes==3
                    trialTypes={'NSS','SSCS','SSNCS'};
                elseif  numTrialTypes==4
                    trialTypes={'LmcNSS','LmncNSS','SSCS','SSNCS'};
                else
                    trialTypes={'NSS','SSCS'};
                end
            elseif strcmp(taskType,'st_saccades')
                trialTypes={'CorrectSTSaccade','FailedSTSaccades'};
            end
            
            %% ECodes 
            % Ecodes Countermanding                         Column
            % 624y / 427y     Fixation cue                  4
            % 644y / 447y     Eye in window                 5
            % 664y / 467y     Fixation cue turned off       6
            % 684y / 487y     Target onset                  7
            % 704y / 507y     Saccade onset or stop signal  8 / 8 or 9, if ecodeout(8)==1503
            % 1030            Reward                        11 / 11 (10 if ecodeout(8)~=1503 ?)
            % 17385           Bad or aborted trial
            % 16386           error code for early saccade  10
            
            % ECodes  Self-timed                                     Column
            %         1001		ENABLECD                                    1
            %         602y		Basecode                                    2
            %         622y		Onset fix target                            4
            %         642y		Eye in window                               5
            %         662y		Flash of cue light                          6
            %         682y		Cue turned off                              7
            %         702y		Rex detected saccade onset                  8
            %         722y		Eye is now in target window                 9
            %         742y		Re-display target after correct trial       10
            %         1030		REWCD                                       11
            %         16386     error code for early saccade                8
            %         17385     Bad or aborted trial                        9

            eventColumns=[8,6,11]; %Saccade/Target/Reward for self timed saccades. [8,7,8,11] for Countermanding
            %             allcodes(trialIdx,:)
            rexnumtrials=size(allcodes,1);
            recDir=[fileDir filesep]; cd(recDir);
            if strcmp(taskType,'gapstop')
                alignedData=GetCmdAlignedData([fileName '_REX.mat'],alignTypes{alignTypeNum},{rexnumtrials,recDir});
            elseif strcmp(taskType,'st_saccades')
                alignedData=GetSTAlignedData([fileName '_REX.mat'],...
                    alignTypes{alignTypeNum},{rexnumtrials,recDir},allcodes);
            end
            [alignedData.rast]=deal(alignedData.rasters);
            [alignedData.alignt]=deal(alignedData.alignidx);
            [alignedData.trialnb]=deal(alignedData.trials);
            [alignedData.evttime]=deal(alignedData.allgreyareas);
            allProcData(fileNum).(alignTypes{alignTypeNum})=rmfield(alignedData,...
                {'dir','trigtosac','sactotrig','trigtovis','vistotrig','eyevel',...
                'sacspecs','stats','alignlabel','savealignname','extras',...
                'rasters','trials','alignidx','allgreyareas'});
            
            %% compare with data from processed file, trial type by trial (rows)
            %         for rowNum=1:numTrialTypes
            %             if (alignTypeNum==1 && rowNum>1) || (alignTypeNum==3 && rowNum~=3)...
            %                     || (alignTypeNum==4 && rowNum==3)
            %                 skipIt=true;
            %             else
            %                 skipIt=false;
            %             end
            %             if ~skipIt
            %                 trialIdx = structData{1,alignTypeNum}(rowNum).trialnb; %structData.(alignTypes{alignTypeNum})(rowNum).trialnb;
            %                 evtCol=eventColumns(alignTypeNum);
            %                 if mode(allcodes(trialIdx,8))==1503 && alignTypeNum==3
            %                     evtCol=evtCol+1;
            %                 end
            %                 if rowNum>1 && mode(allcodes(trialIdx,8))~=1503 && alignTypeNum==4
            %                     mode(allcodes(trialIdx,10))% check if it's 10 instead of 11
            %                 end
            %                 % get spikes
            %                 trialSpikeTimes = allspk(trialIdx,:);
            %                 trialSpikeTimes_alignedToEvent = trialSpikeTimes-alltimes(trialIdx,evtCol);
            %                 minTime=abs(min(min(trialSpikeTimes_alignedToEvent)));
            %                 rasters_fromProc=zeros(size(trialSpikeTimes_alignedToEvent,1),...
            %                     minTime+max(max(trialSpikeTimes_alignedToEvent))+1);
            %
            %                 %% collect data from process data file
            % %                 figure
            %                 ind = 1;
            %                 for i = 1:size(trialSpikeTimes_alignedToEvent,1)
            %                     currentTrial = trialSpikeTimes_alignedToEvent(i,:);
            % %                     for k = 1:length(currentTrial)
            % %                                  plot([currentTrial(k) currentTrial(k)],[ind-1 ind],'k'), hold on
            % %                     end
            %                     rasters_fromProc(i,currentTrial(~isnan(currentTrial))+minTime+1)=1;
            %                     ind = ind+1;
            %                 end
            %
            %                 %                 title('spike times aligned to stop signal (Raw Data)')
            %                 %                 xlim([-2500 400])
            %                 %% get equivalent data from alignData
            % %                                 figure;
            %                 rasters_fromNewProc = alignedData(rowNum).rasters;
            %                 trialAlignment_fromNewProc = alignedData(rowNum).alignidx;
            % %                                 ind = 1;
            % %                                 for i = 1:size(rasters_fromNewProc,1)
            % %                                     currentTrial = rasters_fromNewProc(i,:);
            % %                                     spikesInCurrentTrial = find(currentTrial==1);
            % %                                     spikesInCurrentTrialRelativeToAlignment = spikesInCurrentTrial-trialAlignment_fromNewProc;
            % %                                     spikesForRow = spikesInCurrentTrialRelativeToAlignment; %for brevity I rename variable
            % %
            % %                                     for k = 1:5:length(spikesForRow)
            % %                                         plot([spikesForRow(k) spikesForRow(k)],[ind-1 ind],'k'),hold on
            % %                                     end
            % %                                     ind = ind+1;
            % %                                 end
            %
            %                 %                 xlim([-2500 400])
            %                 %                 title('spike times aligned to stop signal (cmdata)')
            %
            %                 % recData=data.gsdata.allndata{19, 3} %foo.gsdata.allndata{19, 1};
            %                 % start=recData(3).alignt-600;
            %                 % stop=recData(3).alignt+600;
            %                 % crop_rasters = recData(3).rast(:,start:stop);
            %                 % figure
            %                 % [indy, indx] = ind2sub(size(crop_rasters),find(crop_rasters)); %find row and column coordinates of spikes
            %                 % plot([indx';indx'],[indy';indy'+1],'color','k',...
            %                 %     'LineStyle','-','LineWidth',1.8); % plot rasters
            %                 % set(gca,'xlim',[1 length(start:stop)]);
            %
            %
            %                 %% plots
            %                 if alignTypeNum==1
            %                     start=-600; stop=200;
            %                 elseif alignTypeNum==2
            %                     start=-200; stop=400;
            %                 elseif alignTypeNum==3
            %                     start=-400; stop=400;
            %                 elseif alignTypeNum==4
            %                     start=-600; stop=0;
            %                 end
            %                 convsdf_fromstruct=conv_raster(rasters_fromNewProc,30,0,...
            %                     trialAlignment_fromNewProc+start,trialAlignment_fromNewProc+stop);
            %                 convsdf_fromproc=conv_raster(rasters_fromProc,30,0,minTime+start,minTime+stop);
            %                 figure; hold on;
            %                 plot(convsdf_fromstruct)
            %                 plot(convsdf_fromproc)
            %                 title({['File ' fileName];...
            %                     [trialTypes{rowNum} ' trials aligned to ' alignTypes{alignTypeNum}]})
            %                 close all
            %             end
            %         end
        end
        clearvars -except fixcode fixoffcode tgtcode tgtoffcode saccode ...
            stopcode rewcode errcode1 errcode2 errcode3 basecode ...
            dectgtcode decsaccode alignTypes procDir numFiles stdata cmdata ...
            taskType allProcData
        clear global
    end
    save('cDN_stdata_alldata','allProcData','-v7.3'); %top_cortex_cmdata_alldata