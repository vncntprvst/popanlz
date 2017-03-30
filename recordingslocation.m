function recLoc=recordingslocation(fileList)
% compare db recloc with xls file info

userinfo=SetUserDir; 
cd(userinfo.syncdir);

[~, xlsData] = xlsread('Countermanding Data.xlsx');

dbConn = connect2DB('vp_sldata');
dbData=fetch(dbConn,['SELECT e_file,recloc FROM recordings r WHERE r.e_file IN (' ...
    sprintf('''%s'',' ,fileList{1:end-1}) '''' (fileList{end}) ''')']);

identicalLocationValue=false(size(dbData,1),1);
for dbEntry=1:size(dbData,1)
    xlsFileEntry=~cellfun('isempty',strfind(xlsData(:,1),dbData{dbEntry,1}));
    identicalLocationValue(dbEntry)=strcmp(xlsData{xlsFileEntry,2},dbData{dbEntry,2});
end
% find(identicalLocationValue~=true);
recLoc=[dbData mat2cell(identicalLocationValue==true,ones(size(dbData,1),1))];
