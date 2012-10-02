function [sacnr,sacnrsessionlist,sacnrlocationlist,sacnrdepthlist,sacnrtasklist,compart,statscompile]=...
    rec_class(monknum)

if strcmp(getenv('username'),'SommerVD') || strcmp(getenv('username'),'DangerZone')
    directory = 'C:\Data\Recordings\';
else
    directory = 'B:\data\Recordings\';
end
slash = '\';
%monknum=1;

% get number rows
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

%import data
[~,pfilelist] = xlsread([directory 'procdata.xlsx'],monknum,['A2:A' num2str(numrows)]);
sessionlist = xlsread([directory 'procdata.xlsx'],monknum,['B2:B' num2str(numrows)]);
[~,locationlist] = xlsread([directory 'procdata.xlsx'],monknum,['C2:C' num2str(numrows)]);
depthlist = xlsread([directory 'procdata.xlsx'],monknum,['D2:D' num2str(numrows)]);
[~,tasklist] = xlsread([directory 'procdata.xlsx'],monknum,['F2:F' num2str(numrows)]);
activlist = xlsread([directory 'procdata.xlsx'],monknum,['J2:J' num2str(numrows)]);
recnumlist = xlsread([directory 'procdata.xlsx'],monknum,['E2:E' num2str(numrows)]);

%parse data
activlist(isnan(activlist))=0;
%take only analyzed processed files
pfilelist=pfilelist(1:length(activlist));
sessionlist=sessionlist(1:length(activlist));
locationlist=locationlist(1:length(activlist));
depthlist=depthlist(1:length(activlist));
tasklist=tasklist(1:length(activlist));
recnumlist=recnumlist(1:length(activlist));

%narrow down files found to be active
%foo=regexpi(pfilelist,'^\w\d+\w\d\w\d_','match');
or_pfilelist=pfilelist;
pfilelist=cellfun(@(x) x(1:end-1), pfilelist, 'UniformOutput', false);
activefile=or_pfilelist(logical(activlist));
activetfile=pfilelist(logical(activlist));
acsessionlist=sessionlist(logical(activlist));
aclocationlist=locationlist(logical(activlist));
acdepthlist=depthlist(logical(activlist));
actasklist=tasklist(logical(activlist));
maxrecnum=max(recnumlist(logical(activlist)));

%fin out unique neurons (because multiple recordings for each neuron)
[~,sacnridx1]=unique(activetfile);
sacnr= activefile(sacnridx1);
sacnrsessionlist=acsessionlist(sacnridx1);
sacnrlocationlist=aclocationlist(sacnridx1);
sacnrdepthlist=acdepthlist(sacnridx1);
sacnrtasklist=actasklist(sacnridx1);

% crude segmentation: set depth limits for top cortex / dentate / bottom
% cortex
if monknum==1
    cdn_depth=19000;
    bcx_depth=22000;
elseif monknum==2
    cdn_depth=11000;
    bcx_depth=17000;
end
compart=cell(length(sacnrdepthlist),1);
compart(sacnrdepthlist<cdn_depth)={'top_cortex'};
compart(sacnrdepthlist>=cdn_depth & sacnrdepthlist<=bcx_depth)={'dentate'};
compart(sacnrdepthlist>bcx_depth)={'bottom_cortex'};

% get type activity found
algdir=[directory,'processed',slash,'aligned',slash];
algdirlisting=dir(algdir);
algfiles = {algdirlisting(:).name};
algfiles = algfiles(cellfun('isempty',strfind(algfiles,'2SH')));
algfiles = algfiles(~cellfun('isempty',strfind(algfiles,'_sac')));
algfiles = cellfun(@(x) x(1:end-8), algfiles, 'UniformOutput', false);
%making sure all files have their sac alignement file
if length(ismember(sacnr,algfiles'))~=size(sacnr,1)
    disp('missing alignment file')
    return
end

statscompile = struct('p',{},'h',{},'bestmeandiff',{},'recrasters',{});
for getstat=1:size(sacnr,1)
    load([algdir,sacnr{getstat},'_sac.mat']);
    filehstat=arrayfun(@(x) sum(x{:}.h), {dataaligned(~cellfun(@isempty, {dataaligned.stats})).stats});
    if logical(sum(filehstat)) %hopefully it is!
        gooddir=find(filehstat);
        if length(filehstat)<size(dataaligned,2)
            dataaligned=dataaligned(~cellfun(@isempty, {dataaligned.stats}));
        end
        %         pvals=cell(length(find(filehstat)),1);
        %         hvals=cell(length(find(filehstat)),1);
        for getph=1:length(find(filehstat))
            statscompile(getstat,1).p{getph}={dataaligned(gooddir(getph)).stats.p};
            statscompile(getstat,1).h{getph}={dataaligned(gooddir(getph)).stats.h};
            statscompile(getstat,1).recrasters{getph}=dataaligned(1,gooddir(getph)).rasters;
            statscompile(getstat,1).bestmeandiff{getph}=...
                [find(statscompile(getstat,1).p{getph}{:}==...
                max(statscompile(getstat,1).p{getph}{:}(2*(find(statscompile(getstat,1).h{getph}{:})))))/2,...
                max(statscompile(getstat,1).p{getph}{:}(2*(find(statscompile(getstat,1).h{getph}{:}))))];
        end
        
        %         if we want to plot:
        %         figure;
        %         start = dataaligned(1,goodh).alignidx - 500;
        %         stop = dataaligned(1,goodh).alignidx + 300;
        %         hold on;
        %         for trial=1:size(recrasters,1)
        %             spiketimes=find(recrasters(trial,start:stop)); %converting from a matrix representation to a time collection, within selected time range
        %             plot([spiketimes;spiketimes],[ones(size(spiketimes))*trial;ones(size(spiketimes))*trial-1],'k-');
        %         end
        %         sumall=sum(recrasters(:,start:stop));
        %         sdf=spike_density(sumall,12)./size(recrasters,1); %instead of number of trials
        %         set(gca,'TickDir','out'); % draw the tick marks on the outside
        %         set(gca,'YTick', []); % don't draw y-axis ticks
        %         set(gca,'YDir','reverse');
        %         set(gca,'YColor',get(gcf,'Color')); % hide the y axis
        %         set(gca,'XColor',get(gcf,'Color')); % hide the y axis
        %         box off
        %         xlim=get(gca,'xlim');
        %         axes
        %         plot(sdf,'Color','b','LineWidth',2);
        %         set(gca,'Color','none','YAxisLocation','right','TickDir','out', ...
        %             'FontSize',8,'xlim',xlim);
        %         hold off
    end
    for otherfiles=1:maxrecnum-1
        nextfile=[sacnr{getstat}(1:end-1),num2str(str2num(sacnr{getstat}(end-1:end))+otherfiles)];
        if exist([directory,'processed',slash,'aligned',slash,nextfile,'_sac.mat'],'file')
            load([algdir,nextfile,'_sac.mat']);
            filehstat=arrayfun(@(x) sum(x{:}.h), {dataaligned(~cellfun(@isempty, {dataaligned.stats})).stats});
            if logical(sum(filehstat)) %maybe not!
                gooddir=find(filehstat);
                if length(filehstat)<size(dataaligned,2)
                    dataaligned=dataaligned(~cellfun(@isempty, {dataaligned.stats}));
                end
                for getph=1:length(find(filehstat))
                    statscompile(getstat,otherfiles+1).p{getph}={dataaligned(gooddir(getph)).stats.p};
                    statscompile(getstat,otherfiles+1).h{getph}={dataaligned(gooddir(getph)).stats.h};
                    statscompile(getstat,otherfiles+1).recrasters{getph}=dataaligned(1,gooddir(getph)).rasters;
                    statscompile(getstat,otherfiles+1).bestmeandiff{getph}=...
                        [find(statscompile(getstat,otherfiles+1).p{getph}{:}==...
                        max(statscompile(getstat,otherfiles+1).p{getph}{:}(2*(find(statscompile(getstat,otherfiles+1).h{getph}{:})))))/2,...
                        max(statscompile(getstat,otherfiles+1).p{getph}{:}(2*(find(statscompile(getstat,otherfiles+1).h{getph}{:}))))];
                end
            end
        end
    end
end
end


% [~,uniquefilelist] = xlsread([directory 'procdata.xlsx'],4,'A2:A64');
% foo=ismember(uniquefilelist,activefile);
% activeuniquefiles=uniquefilelist(foo);
% look at unique files that were discarded

