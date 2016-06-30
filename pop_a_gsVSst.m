function pop_a_gsVSst(gsdata,stdata,conn)
global directory slash;
if isempty(directory)
    [directory,slash]=SetUserDir;
end

%number cells
gs.goodrecs=~cellfun('isempty',gsdata.allsacdelay);
st.goodrecs=~cellfun('isempty',stdata.allsacdelay);

%remove bad apples
fn = fieldnames(gsdata);
for lp=1:length(fn)
    gsdata.(fn{lp})=gsdata.(fn{lp})(gs.goodrecs,:);
end

gs.cellnum=sum(~cellfun('isempty',gsdata.allsacdelay));
st.cellnum=sum(~cellfun('isempty',stdata.allsacdelay));

disp([num2str(gs.cellnum) ' cells for gs, '  num2str(st.cellnum) ' cells for st'])

queries{1} = ['SELECT a_file FROM sorts s INNER JOIN recordings r on s.recording_id_fk = r.recording_id WHERE sort_id IN (' ...
    sprintf('%.0f,' ,cellfun(@(x) x.sort_id,gsdata.alldb(1:end-1,1))) num2str(gsdata.alldb{end,1}.sort_id) ')'];
% gs.recnames=fetch(conn,query);

queries{2} = ['SELECT a_file FROM sorts s INNER JOIN recordings r on s.recording_id_fk = r.recording_id WHERE sort_id IN (' ...
    sprintf('%.0f,' ,cellfun(@(x) x.sort_id,stdata.alldb(1:end-1,1))) num2str(stdata.alldb{end,1}.sort_id) ')'];
% st.recnames=fetch(conn,query);

[commoncells,fileidx]=compare_db_filelists(queries,conn);

gs.commoncells=find(fileidx{1});
%% individual sessions
% figure; hold on;
% for ccell=1:length(gs.commoncells)
%     sacdelays=gsdata.allsacdelay{gs.commoncells(ccell), 1}{1, 1};
%     xlim=ceil(max(sacdelays)/100)*100;
%     [saclatquant,saclatxlims]=hist(sacdelays,0:25:xlim);
%     saclatfreq=saclatquant./sum(saclatquant);
%     
%     plot(saclatxlims,saclatfreq);
% end
%% all sessions together
allsacdelays=cellfun(@(x) x.nsst,gsdata.allsacdelay,'UniformOutput',false);
%% using only common cells
% allsacdelays=cellfun(@(x) x.nsst,gsdata.allsacdelay(gs.commoncells),'UniformOutput',false);
allsacdelays=[allsacdelays{:}];
xlim=ceil(max(allsacdelays)/100)*100;
[saclatquant,saclatxlims]=hist(allsacdelays,0:25:xlim);
saclatfreq=saclatquant./sum(saclatquant);
figure;plot(saclatxlims,gauss_filtconv(saclatfreq,1));

%% title('Saccade Latency','FontName','calibri','FontSize',15);
hxlabel=xlabel(gca,'Saccade latency','FontName','calibri','FontSize',12);
set(gca,'XTick',0:200:max(get(gca,'xlim')),'TickDir','out','box','off');
hylabel=ylabel(gca,'Proportion of saccades','FontName','calibri','FontSize',12);
curylim=get(gca,'YLim');
set(gca,'TickDir','out','box','off');

%%
st.commoncells=find(fileidx{2}); %st.commoncells=1:size(fileidx{2},1);
properdelay=false(size(st.commoncells,1),1);
% figure; hold on;
for ccell=1:length(st.commoncells)
    sacdelays=stdata.allsacdelay{gs.commoncells(ccell), 1}{1, 1};
% sacdelays=stdata.allsacdelay{ccell}{1, 1};
    if median(sacdelays)<750
%         commoncells{ccell}
        continue
    end
%     xlim=ceil(max(sacdelays)/100)*100;
%     [saclatquant,saclatxlims]=hist(sacdelays,0:25:xlim);
%     saclatfreq=saclatquant./sum(saclatquant);
%     
%     plot(saclatxlims,saclatfreq);
    properdelay(ccell)=true;
end
%% (almost) all sessions together
allsacdelays=cellfun(@(x) x{:},stdata.allsacdelay(properdelay),'UniformOutput',false);
%% only common recordings
% allsacdelays=cellfun(@(x) x{:},stdata.allsacdelay(st.commoncells(properdelay)),'UniformOutput',false);
allsacdelays=[allsacdelays{:}];
allsacdelays=allsacdelays(allsacdelays>800);
xlim=ceil(max(allsacdelays)/100)*100;
[saclatquant,saclatxlims]=hist(allsacdelays,0:25:xlim);
saclatfreq=saclatquant./sum(saclatquant);
figure;plot(saclatxlims,gauss_filtconv(saclatfreq,2));
% foo=gauss_filtconv(saclatfreq,1);figure;plot(saclatxlims,foo)
% foo=fullgauss_filtconv(saclatfreq,2,1);figure;plot(saclatxlims(7:end-6),foo)
%%
title('Saccade Latency','FontName','calibri','FontSize',15);
hxlabel=xlabel(gca,'Saccade latency','FontName','calibri','FontSize',12);
set(gca,'XTick',0:200:max(get(gca,'xlim')),'TickDir','out','box','off');
hylabel=ylabel(gca,'Proportion of saccades','FontName','calibri','FontSize',12);
curylim=get(gca,'YLim');
set(gca,'TickDir','out','box','off');
end