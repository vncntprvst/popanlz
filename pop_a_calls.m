function pop_a_calls(optionNb)

global directory slash;

%% settings
userinfo=UserDirInfo;  %SetUserDir;
directory=userinfo.directory;
slash=filesep;
% userinfo.user,userinfo.dbldir,userinfo.mapdr,userinfo.servrep,userinfo.mapddataf

cd(userinfo.syncdir);
%% load data (preprocessing from pop_analysis_preproc. See alignment parameters line 372. See also prealign > alignrasters)
load('cDn_cmdata.mat'); %cDn_cmdata.mat  top_cortex_cmdata.mat
unitsDBinfo=cmdata.alldb;
CCNdb = connect2DB('vp_sldata');

% get cluster index from database
unitList=unitsDBinfo.unit_id;  %cellfun(@(x) x.unit_id, unitsDBinfo); %changed to a table
recordingList=unitsDBinfo.rec_id; % cellfun(@(x) x.rec_id, unitsDBinfo);
clusterIdx=fetch(CCNdb,['SELECT profile_type FROM clusters c WHERE c.cluster_id IN (' ...
    sprintf('%.0f,' ,unitList(1:end-1)) num2str(unitList(end)) ')']);
fileNames=fetch(CCNdb,['SELECT a_file FROM recordings r WHERE r.recording_id IN (' ...
    sprintf('%.0f,' ,recordingList(1:end-1)) num2str(recordingList(end)) ')']);
[~,resort]=sort(unitList);[~,resort]=sort(resort);clusterIdx=[clusterIdx{resort}]'; 
[~,resort]=sort(recordingList);[~,resort]=sort(resort);cmdata.fileNames=fileNames(resort);

% check behavioral data: keep only ones with correct behavior data
badapl=isnan(clusterIdx) | cellfun('isempty',cmdata.allsacdelay);

fn = fieldnames(cmdata);
for lp=1:length(fn)
    cmdata.(fn{lp})=cmdata.(fn{lp})(~badapl,:);
end

%% Call task specific analysis
%make separate calls for different conditions
calloptions={'compare_st','behavior_values','task_related','multi_task_normalization',...
    'baseline_thd','trial_by_trial','singlessd','allssd_basic_dft',...
    'allssd_basic_multidct','allssd_control'};
call=calloptions{optionNb};

switch call
    case 'behavior_values'
        %simply plot psychometric and tachometric curves
        CmdBehaviorResults(cmdata,CCNdb);
    case 'compare_st'
        load('cDn_stdata.mat');
        pop_a_gsVSst(cmdata,stdata,CCNdb);
    case 'task_related'
        recloc='cDn';
        task='stdata';
        load('cDn_stdata.mat');
        stdata.taskrelated=pop_a_task_related(stdata);
%         %% Save field to data
%         cd(userinfo.syncdir)
%         save([recloc '_' task],task,'-v7.3');
    case 'multi_task_normalization'
        recloc='cDn';
        task='cmdata';
%         load('cDn_stdata.mat');

%number cells
goodrecs=~cellfun('isempty',cmdata.allsacdelay);
% st.goodrecs=~cellfun('isempty',stdata.allsacdelay);

%remove bad apples
% fn = fieldnames(data);
% for lp=1:length(fn)
%     data.(fn{lp})=data.(fn{lp})(goodrecs,:);
% end
        epochs={[0 250];[300 300]};
        cmdata.normData=cell(size(goodrecs,1),5);
        cmdata.normFactor=NaN(size(goodrecs,1),size(epochs,1)+1);
        [cmdata.normData(goodrecs,:),cmdata.normFactor(goodrecs,:)]=pop_a_normalization(cmdata.allndata(goodrecs,:),epochs);
%       %% Save field to data
        cd(userinfo.syncdir)
        save([recloc '_' task],task,'-v7.3');
    case 'baseline_thd'
        recloc='cDn';
        task='cmdata';
%         load('cDn_stdata.mat');
        cmdata.bslThd=pop_a_bsl_threshold(cmdata,CCNdb);
%         %% Save field to data
%         cd(userinfo.syncdir)
%         save([recloc '_' task],task,'-v7.3');
    case 'trial_by_trial'
        [behav, neur]=trialbytrial(cmdata,CCNdb);
    case 'singlessd'
        %% single ssd
        proc_option.recluster=0;    %need to re-cluster?
        
        proc_option.prefdironly=0;  %only prefered direction?
        proc_option.singlessd=1;    %only best SSD?
        proc_option.ssdpkalign=0;   %classify according to when peak occurs with respect to stop signal
        proc_option.popplots=1;     %plot population plots?
        proc_option.basicplots=1;   %plot basic sac/tgt plots?
        proc_option.controlplots=0; %plot control plots (reactive sac/reward)?
        proc_option.defaultplot=1;  %default or multi-condition plots?
        
        proc_option.printplots=0;   %printplots?
        
        pop_a_countermanding(cmdata,proc_option,CCNdb);
    case 'allssd_basic_dft'
        %% all ssd / basic plots / default
        proc_option.recluster=0;    %need to re-cluster?
        
        proc_option.prefdironly=0;  %only prefered direction?
        proc_option.singlessd=0;    %only best SSD?
        proc_option.ssdpkalign=1;   %classify according to when peak occurs with respect to stop signal        
        proc_option.popplots=1;     %plot population plots?
        proc_option.basicplots=0;   %plot basic sac/tgt plots?
        proc_option.controlplots=0; %plot control plots (reactive sac/reward)?
        proc_option.defaultplot=1;  %default or multi-condition plots?
        
        proc_option.printplots=0;   %printplots?
        
%         [compgssdf,clusgsndata,startstop]=Figure4_PopSDFs(cmdata,proc_option,CCNdb);
%         
%         % calculate normalization factor on the spot        
%         dataForNorm= [cellfun(@(x) [x.NSStrial],...
% {compgssdf.alignTgt},'UniformOutput',false)', cellfun(@(x) [x.NSStrial],...
% {compgssdf.alignSac},'UniformOutput',false)'];
%         normFactor=FindNormFactor(dataForNorm, [startstop(2,1) startstop(1,1)]);
%         
%         Figure4_PlotPopSDFs(compgssdf,clusgsndata,proc_option)
%         
        
%         Figure4_PopSDFs_SaccadeAlignment(cmdata,proc_option,CCNdb)
%         Figure4_PopSDFs_StopSignalAlignment(cmdata,proc_option,CCNdb)
        pop_a_countermanding(cmdata,proc_option,CCNdb);
    case 'allssd_basic_multicdt'
        %% all ssd / basic plots / tri conditions
        proc_option.recluster=0;    %need to re-cluster?
        
        proc_option.prefdironly=0;  %only prefered direction?
        proc_option.singlessd=0;    %only best SSD?
        proc_option.ssdpkalign=0;   %classify according to when peak occurs with respect to stop signal
        proc_option.popplots=1;     %plot population plots?
        proc_option.basicplots=1;   %plot basic sac/tgt plots?
        proc_option.controlplots=0; %plot control plots (reactive sac/reward)?
        proc_option.defaultplot=0;  %default or multi-condition plots?
        
        proc_option.printplots=0;   %printplots?
        
        pop_a_countermanding(cmdata,proc_option,CCNdb);
    case 'allssd_control'
        %% all ssd / control plots
        proc_option.recluster=0;    %need to re-cluster?
        
        proc_option.prefdironly=0;  %only prefered direction?
        proc_option.singlessd=0;    %only best SSD?
        proc_option.ssdpkalign=0;   %classify according to when peak occurs with respect to stop signal        
        proc_option.popplots=1;     %plot population plots?
        proc_option.basicplots=0;   %plot basic sac/tgt plots?
        proc_option.controlplots=1; %plot control plots (reactive sac/reward)?
        proc_option.defaultplot=0;  %default or multi-condition plots?
        
        proc_option.printplots=0;   %printplots?
        
        pop_a_countermanding(cmdata,proc_option,CCNdb);
end
end