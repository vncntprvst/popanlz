% similar to rec_class, but after run rec_class and sortstats

function [sacnr,sacnrsessionlist,sacnrlocationlist,sacnrdepthlist,sacnrtasklist,compart,statscompile]=...
    rec_class(monknum)

if strcmp(getenv('username'),'SommerVD') || strcmp(getenv('username'),'DangerZone')
    directory = 'C:\Data\Recordings\';
else
    directory = 'B:\data\Recordings\';
end
slash = '\';


%% Rigel and Sixx data
rigelcxneuron={'R148L7A0_17503'
'R142L7A1_16850'
'R123L5P5_18400'
'R113L6A2_18900'
'R114L6A2_27810'
'R143L7A2_25651'
'R106L7A2_18600'
'R119L7A2_18600'
'R146L6A0_17500'
'R99L7A0_17400'
'R78L6A0_15240'
'R115L6A0_16100'
'R105L7A2_25160'
'R116L6A1_14860'
'R115L6A0_16300'
'R147L6P1_16200'
'R147L6P1_15500'
'R127L5P3_25110'
'R106L7A2_18101'
'R107L6A3_18510'
'R109L6A1_15600'
'R147L6P1_17951'
'R137L6P3_14200'
'R139L6A2_18901'
'R109L6A1_16400'
'R111L6A2_12401'
'R144L7A2_17052'
'R107L6A3_18810'
'R141L6A1_22020'
'R94L6A3_18100'
'R135L6P3_22390'
'R100L7A0_22800'
'R103L7A1_23800'
'R122L5P5_15852'
'R92L6A0_23000'
'R97L7A1_17902'
'R140L6A1_26120'
'R146L6A0_16153'
'R146L6A0_16753'
'R111L6A2_23110'
'R118L7A1_22590'
'R98L7A0_17200'
'R100L7A0_23600'
'R102L7A1_17500'
'R113L6A2_25250'
'R114L6A2_14100'
'R128L5P4_16352'
'R129L5P3_27030'
'R135L6P3_16950'
'R106L7A2_17650'
'R108L6A3_26870'
'R117L7A2_14500'
'R120L7A2_25010'
'R140L6A1_25760'
'R142L7A1_15750'
'R93L6A1_18000'
'R103L7A1_24000'
'R111L6A2_24250'
'R116L6A1_15300'
'R127L5P3_22010'};

sixxcxneuron={'S77L4A2_9751'
'S77L4A2_9750'
'S78L4A3_9651'
'S79L3A3_9502'
'S80L2A3_9901'
'S80L2A3_9102'
'S82L1A8_23301'
'S85L3A5_10412'
'S98L2A5_10705'
'S99L2A5_10301'
'S107L4A4_9151'};

%% 
compallstrinfo=[];
compallnuminfo=[];
for monknum=1:2
%% Get number rows
exl = actxserver('excel.application');
exlWkbk = exl.Workbooks;
exlFile = exlWkbk.Open([directory 'procdata.xlsx']);
exlSheet = exlFile.Sheets.Item(monknum);% e.g.: 2 = Sixx
robj = exlSheet.Columns.End(4);
numrows = robj.row;
% if numrows==1048576 %empty document
%     numrows=1;
% end
Quit(exl);

%% import data
[allnuminfo,allstrinfo] = xlsread([directory 'procdata.xlsx'],monknum,['A2:O' num2str(numrows)]);
allnuminfo=allnuminfo(:,[1 3 4 9 11]);
allstrinfo=allstrinfo(:,[1 3 6 7 8 9 11 13 14 15]);

%% parse data
allnuminfo(isnan(allnuminfo(:,4)),4)=0;

%% narrow down files found to be active
% or_pfilelist=pfilelist;
% pfilelist=cellfun(@(x) x(1:end-1), pfilelist, 'UniformOutput', false);

%allstrinfo=allstrinfo(allnuminfo(:,4)>0,:);
%allnuminfo=allnuminfo(allnuminfo(:,4)>0,:);

%fin out unique neurons (because multiple recordings for each neuron)
%[~,sacnridx1]=unique(activetfile);

%% crude segmentation: set depth limits for top cortex / dentate / bottom
% cortex
% tcxidx=~cellfun(@isempty,regexp(allstrinfo(:,4),'top_cortex'));
% cdnidx=~cellfun(@isempty,regexp(allstrinfo(:,4),'dentate'));
% bcxidx=~cellfun(@isempty,regexp(allstrinfo(:,4),'bottom_cortex'));

%% compile info
compallstrinfo=[compallstrinfo;allstrinfo];
compallnuminfo=[compallnuminfo;allnuminfo];

end
%rigel and sixx cx neuron
% rigelcxneuronidx=ismember(allstrinfo(:,1),rigelcxneuron);
% sixxcxneuronidx=ismember(allstrinfo(:,1),sixxcxneuron);

    rneuronidx=ismember(compallstrinfo(:,1),rigelcxneuron);
    sneuronidx=ismember(compallstrinfo(:,1),sixxcxneuron);
    allneuridx=(rneuronidx | sneuronidx);
    
cxfilelist=compallstrinfo(allneuridx,1);%(tcxidx | bcxidx)
infolist=cell(size(cxfilelist,1),10);
infolist(:,1:8)=compallstrinfo(allneuridx,3:10);%(tcxidx | bcxidx)
infolist(:,9:10)=num2cell(compallnuminfo(allneuridx,4:5));%(tcxidx | bcxidx)

%getting unique neuron list
alluniquen=cellfun(@(x) x(1:end-1), compallstrinfo(:,1), 'UniformOutput', false);
[alluniquen,alluniquenidx]=unique(alluniquen);
%and same for cxneuronlist
allunqcxn=cellfun(@(x) x(1:end-1), cxfilelist, 'UniformOutput', false);
[allunqcxn,allunqcxnidx]=unique(allunqcxn); % we get 70 unique neurons
cxfilelist=cxfilelist(allunqcxnidx);

% make sure all alignment files exist
algdir=[directory,'processed',slash,'aligned',slash];
algdirlisting=dir(algdir);
algfiles = {algdirlisting(:).name};
algfiles = algfiles(cellfun('isempty',strfind(algfiles,'2SH')));
algfiles = algfiles(~cellfun('isempty',strfind(algfiles,'_sac')));
algfiles = cellfun(@(x) x(1:end-8), algfiles, 'UniformOutput', false);
%making sure all files have their sac alignement file
if length(ismember(cxfilelist,algfiles'))~=size(cxfilelist,1)
    disp('missing alignment file')
end
%referencing neurons which have multiple alignment files
ucxfn_task=cell(length(cxfilelist),1);
    ualgfiles=cellfun(@(x) x(1:end-1), algfiles, 'UniformOutput', false);
for ucxfn=1:length(cxfilelist)
    ucxfn_task(ucxfn)={algfiles(ismember(ualgfiles,allunqcxn(ucxfn)))};
end


%% just getting the leadtime, should have done this before
leadtimes=cell(length(cxfilelist),3);
for cxf=1:length(cxfilelist)
    %cxf
    modulationinfo_task=cell(length(ucxfn_task{cxf}),1);
    originalbestfile=find( ismember(ucxfn_task{cxf},cxfilelist));
    for cxf_task=1:length(ucxfn_task{cxf})
        %cxf_task
        load([algdir,ucxfn_task{cxf}{cxf_task},'_sac.mat']);        
        %alldirs=find(~cellfun(@isempty, {dataaligned.alignidx}));
        gooddirs=find(arrayfun(@(x) nansum(x{:}.h), {dataaligned(~cellfun(@isempty, {dataaligned.stats})).stats}));
        maxmeandiffs=arrayfun(@(x) max(x{:}.p), {dataaligned(gooddirs).stats});
        bestdirs=gooddirs(maxmeandiffs==max(maxmeandiffs));
        if ~isempty(bestdirs)
        modulationinfo=cell(length(gooddirs),3);

        [~,~,~,modulationinfo]=findauc(cxfilelist{cxf},dataaligned,'active');
        
        %modulationinfo is: timetoonset,aligntime-ppkt (bascically, is peak time pre or post-sac) ,slopeshape;
        modulationinfo=[cell(length(gooddirs),2) modulationinfo];
        modulationinfo(bestdirs,1:2)={'best',bestdirs};
        modulationinfo(gooddirs(gooddirs~=bestdirs),1)={'rest'};
        otherdirs=gooddirs(gooddirs~=bestdirs);
        for othd=1:length(otherdirs)
            modulationinfo(otherdirs(othd),2)={otherdirs(othd)};
        end
        modulationinfo_task{cxf_task}=modulationinfo;
        else
            modulationinfo_task{cxf_task}=[];
        end
    end
    % to keep referencing correct, sort excel rows by sessions!
    leadtimes{cxf,1}=modulationinfo_task;
    % even so ucxfn_task is sorted, not compallstrinfo 
    [~,matchfidx]=sort(compallstrinfo(ismember(compallstrinfo(:,1),ucxfn_task{cxf}),1));
    findtasks=compallstrinfo(ismember(compallstrinfo(:,1),ucxfn_task{cxf}),3);
    leadtimes{cxf,2}=findtasks(matchfidx);
    leadtimes{cxf,3}=originalbestfile;
end

%% plotting leadtimes histogram
% find cells with multiple recordings
multitaskcxn=cellfun(@(x) sum(~cellfun('isempty',x)), leadtimes);
multtleadtimes=leadtimes(multitaskcxn(:,2)>1,:);

%separate by tasks
stsacidx=cellfun(@(x) strcmp('st_saccades',x), leadtimes(:,2),'UniformOutput', false);
tokidx=cellfun(@(x) strcmp('tokens',x), leadtimes(:,2),'UniformOutput', false);
gapstopidx=cellfun(@(x) strcmp('gapstop',x), leadtimes(:,2),'UniformOutput', false);
optlocidx=cellfun(@(x) strcmp('optiloc',x), leadtimes(:,2),'UniformOutput', false);

stsaclt=nan(length(stsacidx),1);
for stsaccxn=1:length(stsacidx)  
    if sum(stsacidx{stsaccxn}) % previously: stsacidx{stsaccxn}(min(leadtimes{stsaccxn,3},length(leadtimes{stsaccxn,2})))
        % which was a control to only take leadtimes from task-associated
        % recording which was the supposed best. That didn't work out well.
        ltinfo=leadtimes{stsaccxn}(stsacidx{stsaccxn}); %{min(leadtimes{stsaccxn,3},length(leadtimes{stsaccxn,2}))};
        if ~isempty(ltinfo{:}) 
        bestdirld=cellfun(@(x) strcmp(x(:,1),'best'), ltinfo, 'UniformOutput', false);
        if logical(sum(cell2mat(ltinfo{:}(bestdirld{:},3))))
            stsaclt(stsaccxn)=cell2mat(ltinfo{:}(bestdirld{:},3));
        end
        end
    end
end
stsaclt=stsaclt(stsaclt>-70);

toklt=nan(length(tokidx),1);
for tokcxn=1:length(tokidx)
    %tokcxn
    if sum(tokidx{tokcxn}) % previously: tokidx{tokcxn}(min(leadtimes{tokcxn,3},length(leadtimes{tokcxn,2})))
        % which was a control to only take leadtimes from task-associated
        % recording which was the supposed best. That didn't work out well.
        ltinfo=leadtimes{tokcxn}(tokidx{tokcxn}); %{min(leadtimes{tokcxn,3},length(leadtimes{tokcxn,2}))};
        try 
            iscell(ltinfo{:});
        catch
            continue
        end
        ltinfo=ltinfo(~cellfun('isempty',ltinfo));
        if ~isempty(ltinfo) && logical(sum(ltinfo{1}{3}))
            
            if ~isempty(ltinfo{:})
                bestdirld=cellfun(@(x) strcmp(x(:,1),'best'), ltinfo, 'UniformOutput', false);
                if logical(sum(cell2mat(ltinfo{:}(bestdirld{:},3))))
                    toklt(tokcxn)=cell2mat(ltinfo{:}(bestdirld{:},3));
                end
            end
        end
        end
end
toklt=toklt(toklt>-70);

gapstoplt=nan(length(gapstopidx),1);
for gapstopcxn=1:length(gapstopidx)
    %gapstopcxn
    if sum(gapstopidx{gapstopcxn}) % previously: gapstopidx{gapstopcxn}(min(leadtimes{gapstopcxn,3},length(leadtimes{gapstopcxn,2})))
        % which was a control to only take leadtimes from task-associated
        % recording which was the supposed best. That didn't work out well.
        ltinfo=leadtimes{gapstopcxn}(gapstopidx{gapstopcxn}); %{min(leadtimes{gapstopcxn,3},length(leadtimes{gapstopcxn,2}))};
        try
            iscell(ltinfo{:});
        catch
            continue
        end
            ltinfo=ltinfo(~cellfun('isempty',ltinfo));
            
            if ~isempty(ltinfo) && logical(sum(ltinfo{1}{3}))
                
                if ~isempty(ltinfo{:})
                    bestdirld=cellfun(@(x) strcmp(x(:,1),'best'), ltinfo, 'UniformOutput', false);
                    if logical(sum(cell2mat(ltinfo{:}(bestdirld{:},3))))
                        gapstoplt(gapstopcxn)=cell2mat(ltinfo{:}(bestdirld{:},3));
                    end
                end
        end
    end
end
gapstoplt=gapstoplt(gapstoplt>-70);

optloclt=nan(length(optlocidx),1);
for optloccxn=1:length(optlocidx)
    %optloccxn
    if sum(optlocidx{optloccxn}) % previously: optlocidx{optloccxn}(min(leadtimes{optloccxn,3},length(leadtimes{optloccxn,2})))
        % which was a control to only take leadtimes from task-associated
        % recording which was the supposed best. That didn't work out well.
        ltinfo=leadtimes{optloccxn}(optlocidx{optloccxn}); %{min(leadtimes{optloccxn,3},length(leadtimes{optloccxn,2}))};
        try
            iscell(ltinfo{:});
        catch
            continue
        end
            ltinfo=ltinfo(~cellfun('isempty',ltinfo));
            
            if ~isempty(ltinfo) && logical(sum(ltinfo{1}{3}))
                
                if ~isempty(ltinfo{:})
                    bestdirld=cellfun(@(x) strcmp(x(:,1),'best'), ltinfo, 'UniformOutput', false);
                    if logical(sum(cell2mat(ltinfo{:}(bestdirld{:},3))))
                        optloclt(optloccxn)=cell2mat(ltinfo{:}(bestdirld{:},3));
                    end
                end
        end
    end
end
optloclt=optloclt(optloclt>-70);

% plot overlayed histograms
histfh=figure;
%colormap(jet(4));
xlimvalues=(-round(max([optloclt; gapstoplt; toklt; stsaclt])/100)*100:25:70);
hist(-stsaclt,xlimvalues);
hold on;
hist(-gapstoplt,xlimvalues);
%hist(-toklt,xlimvalues);
hist(-optloclt,xlimvalues);
legendh=legend('Self Timed Saccade Task','Countermanding Task','Optimal Location','Location','NorthWest');
lddatah = findobj(gca,'Type','patch');
% display(lddatah)
%colormap summer
cc=lines(3);
set(lddatah(1),'Facecolor',cc(1,:),'EdgeColor',cc(1,:),'FaceAlpha', 0.1,'LineWidth',1.8);
set(lddatah(2),'Facecolor',cc(2,:),'EdgeColor',cc(2,:),'FaceAlpha', 0.1,'LineWidth',1.8);
set(lddatah(3),'Facecolor',cc(3,:),'EdgeColor',cc(3,:),'FaceAlpha', 0.1,'LineWidth',1.8);
%set(lddatah(4),'Facecolor',cc(4,:),'EdgeColor',cc(4,:),'FaceAlpha', 0.1,'LineWidth',1.8);
        set(gca,'Color','none','TickDir','out','FontSize',8,'FontName','calibri','box','off');
                ylabel(gca,'Number of units','FontName','calibri','FontSize',8);
        xlabel(gca,'Lead Time (ms)','FontName','calibri','FontSize',8);
            title('Distribution of saccadic lead times in lateral cerebellar neurons, by task','FontName','calibri','FontSize',11);
set(legendh,'Interpreter','none', 'Box', 'off','LineWidth',1.5,'FontName','calibri','FontSize',10);




%gather lead time data
% lddata(1,:)=hist(-stsaclt,20);
% lddata(2,:)=hist(-toklt,20);
% lddata(3,:)=hist(-gapstoplt,20);
% lddata(4,:)=hist(-optloclt,20);
% 
% lddatah = bar(lddata','hist');
% set(get(lddatah(1),'BaseLine'),'LineWidth',2,'LineStyle',':')
% %set(gca,'YLim',[1 max(lddatah)]);
% colormap summer % Change the color scheme

%% get best leadtimes
bestlt=cell2mat(cellfun(@(x) x(cellfun(@(x) strcmp('best',x), leadtimes{1}(:,1)),3),leadtimes))
hist(-bestlt)
%% plotting allfiles

effectcat='allcx_tokeep';
SummaryPlot(effectcat,cxfilelist,infolist);

close all force;

