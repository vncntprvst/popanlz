function bslThd=pop_a_bsl_threshold(data,conn)
%get threshold from baseline for each pair of diff activity
% (CSST vs NSST, NCSST vs NSST)
global directory slash;
if isempty(directory)
    [directory,slash]=SetUserDir;
end

%number cells
gs.goodrecs=~cellfun('isempty',data.allsacdelay);
% st.goodrecs=~cellfun('isempty',stdata.allsacdelay);

%remove bad apples
fn = fieldnames(data);
for lp=1:length(fn)
    data.(fn{lp})=data.(fn{lp})(gs.goodrecs,:);
end

% gs.cellnum=sum(~cellfun('isempty',gsdata.allsacdelay));
% st.cellnum=sum(~cellfun('isempty',stdata.allsacdelay));
% 
% disp([num2str(gs.cellnum) ' cells for gs, '  num2str(st.cellnum) ' cells for st'])
% 
% queries{1} = ['SELECT a_file FROM sorts s INNER JOIN recordings r on s.recording_id_fk = r.recording_id WHERE sort_id IN (' ...
%     sprintf('%.0f,' ,cellfun(@(x) x.sort_id,gsdata.alldb(1:end-1,1))) num2str(gsdata.alldb{end,1}.sort_id) ')'];
% % gs.recnames=fetch(conn,query);
% 
% queries{2} = ['SELECT a_file FROM sorts s INNER JOIN recordings r on s.recording_id_fk = r.recording_id WHERE sort_id IN (' ...
%     sprintf('%.0f,' ,cellfun(@(x) x.sort_id,stdata.alldb(1:end-1,1))) num2str(stdata.alldb{end,1}.sort_id) ')'];
% % st.recnames=fetch(conn,query);
% 
% [~,fileidx]=compare_db_filelists(queries,conn);
% 
% gs.commoncells=find(fileidx{1});

%% Convolve traces. 
% SST vs NSS saccade task condition
sigma=15;
bslActivity=cellfun(@(x) {conv_raster(x(1).rast,sigma,x(1).alignt-(600+sigma*3),x(1).alignt+(sigma*3)); ...
    conv_raster(x(2).rast,sigma,x(2).alignt-(600+sigma*3),x(2).alignt+(sigma*3)); ...
    conv_raster(x(3).rast,sigma,x(3).alignt-(600+sigma*3),x(3).alignt+(sigma*3))},...
    data.allndata(:,2), 'UniformOutput',false); % 600ms epoch

% % remove recordings with one trial
% properLengthRecs=cellfun(@(x) length(x)>1,bslActivity);
% gs.goodrecs(gs.goodrecs)=properLengthRecs;
% bslActivity=bslActivity(properLengthRecs);

%% Find baseline threshold 
% mean diff1 / std diff1 / mean diff2 / std diff2
bslThd=nan(size(gs.goodrecs,1),4);
bslThd(gs.goodrecs,1:4)=cell2mat(cellfun(@(bslEpoch) [nanmean(abs(bslEpoch{2}-bslEpoch{1})) ...
    2*nanstd(abs(bslEpoch{2}-bslEpoch{1})) nanmean(abs(bslEpoch{3}-bslEpoch{1})) ...
    2*nanstd(abs(bslEpoch{3}-bslEpoch{1}))],bslActivity, 'UniformOutput',false));
