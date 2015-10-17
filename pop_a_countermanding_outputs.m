
global directory slash;

%% settings
[directory,slash,user,dbldir,mapdr,servrep,mapddataf]=SetUserDir;
CCNdb = connect2DB('vp_sldata');

cd('E:\BoxSync\Box Sync\Home Folder vp35\Sync\SommerLab\projects\countermanding\popclusters\')
load('countermanding_cDn_gsdata.mat');

%% Call task specific analysis
%make separate calls for different conditions
calloptions={'singlessd','allssd_basic_dft','allssd_basic_multidct','allssd_control'};
call=calloptions{2};

switch call
    case 'singlessd'
    %% single ssd
    proc_option.recluster=0;    %need to re-cluster? 

    proc_option.prefdironly=0;  %only prefered direction?
    proc_option.singlessd=1;    %only best SSD?

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

    proc_option.popplots=1;     %plot population plots?
    proc_option.basicplots=1;   %plot basic sac/tgt plots?
    proc_option.controlplots=0; %plot control plots (reactive sac/reward)?
    proc_option.defaultplot=1;  %default or multi-condition plots?

    proc_option.printplots=1;   %printplots?

    pop_a_countermanding(gsdata,proc_option,CCNdb);

    case 'allssd_basic_multicdt'
    %% all ssd / basic plots / tri conditions
    proc_option.recluster=0;    %need to re-cluster? 

    proc_option.prefdironly=0;  %only prefered direction?
    proc_option.singlessd=0;    %only best SSD?

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

    proc_option.popplots=1;     %plot population plots?
    proc_option.basicplots=0;   %plot basic sac/tgt plots?
    proc_option.controlplots=1; %plot control plots (reactive sac/reward)?
    proc_option.defaultplot=0;  %default or multi-condition plots?

    proc_option.printplots=0;   %printplots?

    pop_a_countermanding(gsdata,proc_option,CCNdb);

end