function [ratios,cancelTime,pValues]=precancel_activdiff
% returns ratios,cancelTime

global directory slash;

%% settings
userinfo=SetUserDir;
directory=userinfo.directory;
slash=userinfo.slash;

cd(userinfo.syncdir);
load('cDn_gsdata.mat'); %cDn_gsdata.mat  top_cortex_gsdata.mat

    %remove bad apples
    goodrecs=~cellfun('isempty',gsdata.allsacdelay);

    fieldName = fieldnames(gsdata);
    for lp=1:length(fieldName)
        gsdata.(fieldName{lp})=gsdata.(fieldName{lp})(goodrecs,:);
    end
    
    % get unit cluster info and profiles
    CCNdb = connect2DB('vp_sldata');
    unit_ids=cellfun(@(x) x.unit_id,gsdata.alldb);
    
    [sorted_unit_ids,sunitid_idx]=sort(unit_ids);
    query = ['SELECT c.profile, c.profile_type FROM clusters c WHERE cluster_id IN (' sprintf('%.0f,' ,sorted_unit_ids(1:end-1)') num2str(sorted_unit_ids(end)) ')'];
    profiles = fetch(CCNdb,query);
    sunitid_revidx(sunitid_idx)=1:length(unit_ids);
    clusidx=[profiles{sunitid_revidx,2}];
    clustypes={profiles{sunitid_revidx,1}};
    
    % keep only one cluster
    clus1Idx=clusidx==101;
    for lp=1:length(fieldName)
        gsdata.(fieldName{lp})=gsdata.(fieldName{lp})(clus1Idx,:);
    end

for recNum=1:size(gsdata.allndata,1)

    allSSDs=[gsdata.allssds{recNum, 1}{1, 1};gsdata.allssds{recNum, 1}{1, 2}];
    [~,~,uniqSSDsIdx]=unique(floor(allSSDs/5)*5);
    
    if numel(gsdata.allmssrt_tacho{recNum, 1})>1
        SSRT=round(gsdata.allmssrt_tacho{recNum, 1}{1});
    else
        ratios{recNum,1}=[];
        cancelTime{recNum,1}=[];
        continue;
    end
    for SSDnum=1:length(unique(uniqSSDsIdx))
        try
            prevssd=mode(allSSDs(uniqSSDsIdx==SSDnum));
            
            % CSST Latency matched NSS trials
            % Keeping NSS trials with sac latencies long enough
            % that they would have occured after the stop-signal and ssrt
            alignedata=gsdata.allndata{recNum, 2};
            %keep only first two conditions
            alignedata=alignedata(1:2);
            
            numrast=length(alignedata);
    
            fieldName = fieldnames(alignedata);
            fieldName=fieldName(structfun(@(x) numel(x)>1 & ~ischar(x), alignedata(1, 1),'UniformOutput', true));
            
            for lp=1:length(fieldName)
                try
                    alignedata(1).(fieldName{lp})=alignedata(1).(fieldName{lp})(...
                        gsdata.allsacdelay{recNum, 1}.nsst>prevssd+SSRT,:);
                catch
                    alignedata(1).(fieldName{lp})=alignedata(1).(fieldName{lp})(...
                        :,gsdata.allsacdelay{recNum, 1}.nsst>prevssd+SSRT);
                end
            end
            
            %keep relevant SS trials (with SSD within +-3ms of prevalent SSD)
            %CS trials
            for lp=1:length(fieldName)
                try
                    alignedata(2).(fieldName{lp})=alignedata(2).(fieldName{lp})(logical(arrayfun(@(x) sum(prevssd<=x+3 & prevssd>=x-3),...
                        gsdata.allssds{recNum, 1}{1, 1})),:);
                catch
                    alignedata(2).(fieldName{lp})=alignedata(2).(fieldName{lp})(:,logical(arrayfun(@(x) sum(prevssd<=x+3 & prevssd>=x-3),...
                        gsdata.allssds{recNum, 1}{1, 1})));
                end
            end
            
            %plot related variables
            startCut=400;
            stopCut=600;
            fsigma=15;
            causker=0;
            chrono=1;
            clear sdf;
            for rastnum=1:numrast
                try
                    rasters=alignedata(rastnum).rast;
                    alignidx=alignedata(rastnum).alignt;
                    
                    start=alignidx - startCut;
                    stop=alignidx + stopCut;
                    
                    if start < 1
                        start = 1;
                        disp(['shorten raster ' num2str(recNum)])
                    end
                    if stop > length(rasters)
                        stop = length(rasters);
                        disp(['shorten raster ' num2str(recNum)])
                    end
                    
                    cut_rasters = rasters(:,start:stop); % Isolate rasters of interest
                    isnantrial = isnan(sum(cut_rasters,2)); % Identify nantrials
                    
                    if size(rasters(~isnantrial,:),1)<5 %if less than 5 good trials
                        %useless plotting this
                        sumall=NaN;
                    else
                        sumall=sum(rasters(~isnantrial,start-fsigma*3:stop+fsigma*3));
                    end
                    %     sdf=spike_density(sumall,fsigma)./length(find(~isnantrial)); %instead of number of trials
                    sdf{rastnum}=fullgauss_filtconv(sumall,fsigma,causker)./length(find(~isnantrial)).*1000;
      
                catch
                    sdf{rastnum}=[];
                    continue
                end
            end
%             figure;hold on;
%             plot(sdf{1});
%             plot(sdf{2});

            %% get ratio and cancellation times
            %ratio of activity in 40ms preceding SSRT
            
            if isempty(sdf{1}) || isempty(sdf{2})
                disp('empty sdf');
                ratios{recNum,SSDnum}=[];
                cancelTime{recNum,SSDnum}=[];
                continue
            end
            if length(sdf{1})~=length(sdf{2})
                abs(length(sdf{1})~=length(sdf{2}));
            end
            
            %ratio 
            cancelWindow=startCut+prevssd+SSRT-40:startCut+prevssd+SSRT-1;
            ratios{recNum,SSDnum}=nanmean(sdf{1}(cancelWindow))/nanmean(sdf{2}(cancelWindow));
            if ~isempty(ratios{recNum,SSDnum}) && sum(size(sdf{1})==size(sdf{2}))==2
                [~, ~, pValues{recNum,SSDnum}] = statcond({sdf{1}(cancelWindow) sdf{2}(cancelWindow)},...
                    'method', 'perm', 'naccu', 5000);
            else 
                pValues{recNum,SSDnum}=[];
            end
            % cancelation time
            % differential spike density function exceeded by 2 SD the mean
            % difference in activity during the 600-ms interval before 
            % target presentation, provided the difference reaches 6 SD and
            % remains >2 SD threshold for 50 ms.
            diff_sdf=sdf{1}-sdf{2};
            threshold=nanmean(diff_sdf(1:alignidx-start))+nanstd(diff_sdf(1:alignidx-start)); 
            aboveThdEpochs=bwlabel(diff_sdf(alignidx-start:end)>threshold);
            if max(aboveThdEpochs)>0
                for epoch=1:max(aboveThdEpochs)
                    if logical(sum(diff_sdf(startCut+find(aboveThdEpochs==epoch))>3*threshold)) &&...
                            sum(aboveThdEpochs==epoch)>=50
                        cancelTime{recNum,SSDnum}(epoch)=startCut+prevssd+SSRT-...
                            (startCut+find(aboveThdEpochs==epoch,1));
                    else
                        cancelTime{recNum,SSDnum}(epoch)=nan;
                    end
                end
            end
            
        catch
            ratios{recNum,SSDnum}=[];
            cancelTime{recNum,SSDnum}=[];
        end
    end
end

figure;
histogram([ratios{:}]);
figure;
binEdges=linspace(floor(min([cancelTime{:}])/10)*10,ceil(max([cancelTime{:}])/10)*10,...
    ((ceil(max([cancelTime{:}])/10)*10-floor(min([cancelTime{:}])/10)*10)/10)+1);
histogram([cancelTime{:}],binEdges)
% pValues=[pValues{:}];
% pValues(pValues<0.05)
