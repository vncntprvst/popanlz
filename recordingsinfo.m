function recInfo=recordingsinfo

userinfo=SetUserDir; 
cd(userinfo.syncdir);

datasets={'cDn_gsdata','top_cortex_gsdata'};
[datasetNum,okgo] = listdlg('PromptString','Select dataset to load:',...
    'SelectionMode','single',...
    'ListString',datasets,...
    'OKString','Load');
if okgo
    dataset=datasets{datasetNum};
else %default
    dataset='cDn_gsdata';
end

data=load([dataset '.mat']); %cDn_gsdata.mat  top_cortex_gsdata.mat
dataField=cell2mat(fieldnames(data));
unitsDBinfo=data.(dataField).alldb;
dbConn = connect2DB('vp_sldata');
% get cluster index from database
unitList=cellfun(@(x) x.unit_id, unitsDBinfo);
clusterIdx=fetch(dbConn,['SELECT profile_type, sort_id_fk FROM clusters c WHERE c.cluster_id IN (' ...
    sprintf('%.0f,' ,unitList(1:end-1)) num2str(unitList(end)) ')']);
[~,resort]=sort(unitList);[~,resort]=sort(resort);clusterIdx=cell2mat(clusterIdx(resort,:));
% clusterIdx=[clusterIdx{resort}]';
% check behavioral data: keep only ones with correct behavior data
badapl=isnan(clusterIdx(:,1)) | cellfun('isempty',data.(dataField).allsacdelay);
%remove bad apples
sortsIDs=clusterIdx(~badapl,2);

recIDs = fetch(dbConn, ['SELECT recording_id_fk FROM sorts s WHERE s.sort_id IN (' ...
    sprintf('%.0f,' ,sortsIDs(1:end-1)) num2str(sortsIDs(end)) ')']);
[~,resort]=sort(sortsIDs);[~,resort]=sort(resort);recIDs=cell2mat(recIDs(resort,:));

recInfo=fetch(dbConn,['SELECT e_file,lm_coord,ap_coord,depth,sessions_id_fk,' ...
    'grid_id_fk FROM recordings r WHERE r.recording_id IN (' ...
    sprintf('%.0f,' ,recIDs(1:end-1)) num2str(recIDs(end)) ')']);
[~,resort]=sort(recIDs);[~,resort]=sort(resort);recInfo=recInfo(resort,:);

[sessionIDs,~,allIdx]=unique(cell2mat(recInfo(:,5)));
ugridInfo=fetch(dbConn,['SELECT Subject,Grid_location,Grid_rotated FROM sessions s WHERE s.sessions_id IN (' ...
    sprintf('%.0f,' ,sessionIDs(1:end-1)) num2str(sessionIDs(end)) ')']);
gridInfo=cell(size(recInfo,1),3);
for gridIdx=1:size(recInfo,1)
        gridInfo(gridIdx,:)=ugridInfo(allIdx(gridIdx),:);    
end

% Filename	Monkey	Chamber Location	M-L	A-P	Chamber Rotation   Depth
recInfo=[recInfo(:,1) gridInfo(:,1:2) recInfo(:,2:3) gridInfo(:,3) recInfo(:,4)];

