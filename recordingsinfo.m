function [recCoordInfo,badapl]=recordingsinfo
% from a defined dataset returns three columns:
% file name / recording location / whether location is congruent 
% between database and curated list in xls file

userinfo=SetUserDir; 
cd(userinfo.syncdir);

datasets={'cDn_cmdata','top_cortex_cmdata'};
[datasetNum,okgo] = listdlg('PromptString','Select dataset to load:',...
    'SelectionMode','single',...
    'ListString',datasets,...
    'OKString','Load');
if okgo
    dataset=datasets{datasetNum};
else %default
    dataset='cDn_cmdata';
end

data=load([dataset '.mat']); %cDn_cmdata.mat  top_cortex_cmdata.mat
dataField=cell2mat(fieldnames(data));
unitsDBinfo=data.(dataField).alldb;
dbConn = connect2DB('vp_sldata');
unitList=cellfun(@(x) x.unit_id, unitsDBinfo);
% get cluster index from database
clusterIdx=getclusterindex(dbConn,unitList,'sort_id_fk');

% clusterIdx=[clusterIdx{resort}]';
% check behavioral data: keep only ones with correct behavior data
badapl=isnan(clusterIdx(:,1)) | cellfun('isempty',data.(dataField).allsacdelay);
%remove bad apples
clusterIdx=clusterIdx(~badapl,:);

sortsIDs=clusterIdx(:,2);
recIDs = fetch(dbConn, ['SELECT recording_id_fk FROM sorts s WHERE s.sort_id IN (' ...
    sprintf('%.0f,' ,sortsIDs(1:end-1)) num2str(sortsIDs(end)) ')']);
[~,resort]=sort(sortsIDs);[~,resort]=sort(resort);recIDs=cell2mat(recIDs(resort,:));

recCoordInfo=fetch(dbConn,['SELECT e_file,lm_coord,ap_coord,depth,sessions_id_fk,' ...
    'grid_id_fk FROM recordings r WHERE r.recording_id IN (' ...
    sprintf('%.0f,' ,recIDs(1:end-1)) num2str(recIDs(end)) ')']);
[~,resort]=sort(recIDs);[~,resort]=sort(resort);recCoordInfo=recCoordInfo(resort,:);

[sessionIDs,~,allIdx]=unique(cell2mat(recCoordInfo(:,5)));
ugridInfo=fetch(dbConn,['SELECT Subject,Grid_location,Grid_rotated FROM sessions s WHERE s.sessions_id IN (' ...
    sprintf('%.0f,' ,sessionIDs(1:end-1)) num2str(sessionIDs(end)) ')']);
gridInfo=cell(size(recCoordInfo,1),3);
for gridIdx=1:size(recCoordInfo,1)
        gridInfo(gridIdx,:)=ugridInfo(allIdx(gridIdx),:);    
end

% Filename	Monkey	Chamber Location	M-L	A-P	Chamber Rotation   Depth   ClusterID
recCoordInfo=[recCoordInfo(:,1) gridInfo(:,1:2) recCoordInfo(:,2:3) gridInfo(:,3) recCoordInfo(:,4) mat2cell(clusterIdx(:,1),ones(size(clusterIdx,1),1))];
recCoordInfo=cell2table(recCoordInfo,'VariableNames',{'fileName' 'subject' 'coordinates' 'ML' 'AP' 'rotation' 'depth' 'classification'});

