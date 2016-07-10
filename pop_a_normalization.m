function normFactor=pop_a_normalization(gsdata,stdata,conn)
global directory slash;
if isempty(directory)
    [directory,slash]=SetUserDir;
end

%number cells
gs.goodrecs=~cellfun('isempty',gsdata.allsacdelay);
% st.goodrecs=~cellfun('isempty',stdata.allsacdelay);

%remove bad apples
fn = fieldnames(gsdata);
for lp=1:length(fn)
    gsdata.(fn{lp})=gsdata.(fn{lp})(gs.goodrecs,:);
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
% use rasters from NSS saccade task condition
% motor epoch:  300 ms before to 300ms after saccade onset
% visual epoch: 0–250 ms after target onset

sigma=15;
sacActivity=cellfun(@(x) conv_raster(x(1,1).rast,sigma,x(1,1).alignt-(300+sigma*3),x(1,1).alignt+(299+sigma*3)), gsdata.allndata(:,1), 'UniformOutput',false); %600ms epoch
tgtActivity=cellfun(@(x) conv_raster(x(1,1).rast,sigma,x(1,1).alignt-sigma*3,x(1,1).alignt+(249+sigma*3)), gsdata.allndata(:,2), 'UniformOutput',false); % 250ms epoch

% remove recordings with one trial
properLengthRecs=cellfun(@(x) length(x)>1,sacActivity);
gs.goodrecs(gs.goodrecs)=properLengthRecs;
sacActivity=sacActivity(properLengthRecs);
tgtActivity=tgtActivity(properLengthRecs);

%% Find normalization factor 
% Get peak firing rates during the visual epoch and the motor epoch 
% The larger value is the normalization factor of the cell. 
normFactor=nan(size(gs.goodrecs,1),1);
normFactor(gs.goodrecs)=cellfun(@(visEpoch,motEpoch) max([visEpoch,motEpoch]),tgtActivity,sacActivity, 'UniformOutput',true);

