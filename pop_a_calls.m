function pop_a_calls(optionNb)

global directory slash;

%% settings
userinfo=SetUserDir;
directory=userinfo.directory;
slash=userinfo.slash;
% userinfo.user,userinfo.dbldir,userinfo.mapdr,userinfo.servrep,userinfo.mapddataf
CCNdb = connect2DB('vp_sldata');

cd(userinfo.syncdir);
load('cDn_gsdata.mat'); %cDn_gsdata.mat  top_cortex_gsdata.mat

%% Call task specific analysis
%make separate calls for different conditions
calloptions={'compare_st','trial_by_trial','singlessd','allssd_basic_dft',...
    'allssd_basic_multidct','allssd_control','behavior_values'};
call=calloptions{optionNb};

switch call
    case 'behavior_values'
        %simply plot psychometric and tachometric curves
        CmdBehaviorResults(gsdata,CCNdb);
    case 'compare_st'
        load('cDn_stdata.mat');
        pop_a_gsVSst(gsdata,stdata,CCNdb);
    case 'trial_by_trial'
        [behav, neur]=trialbytrial(gsdata,CCNdb);
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
        
        pop_a_countermanding(gsdata,proc_option,CCNdb);
    case 'allssd_basic_dft'
        %% all ssd / basic plots / default
        proc_option.recluster=0;    %need to re-cluster?
        
        proc_option.prefdironly=0;  %only prefered direction?
        proc_option.singlessd=0;    %only best SSD?
        proc_option.ssdpkalign=0;   %classify according to when peak occurs with respect to stop signal        
        proc_option.popplots=1;     %plot population plots?
        proc_option.basicplots=1;   %plot basic sac/tgt plots?
        proc_option.controlplots=0; %plot control plots (reactive sac/reward)?
        proc_option.defaultplot=1;  %default or multi-condition plots?
        
        proc_option.printplots=0;   %printplots?
        
        pop_a_countermanding(gsdata,proc_option,CCNdb);
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
        
        pop_a_countermanding(gsdata,proc_option,CCNdb);
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
        
        pop_a_countermanding(gsdata,proc_option,CCNdb);
end
end