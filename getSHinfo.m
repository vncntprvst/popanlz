function results=getSHinfo(directory,export)

SHdir=[directory,'SH',slash];
SHdirlisting = dir(SHdir);

%%  extract info

SHdirfileNames = {SHdirlisting.name};

txtfiles=regexpi(SHdirfileNames,'\w*_info.txt$','match');
txtfiles=txtfiles(~cellfun('isempty',txtfiles));

SHfinfo=cellfun(@(x) openreadclose(x,SHdir),txtfiles);

GCLidx=logical(cellfun(@(x) sum(~cellfun('isempty',regexpi(x,'GCL','match'))),SHfinfo)) | ...
    logical(cellfun(@(x) sum(~cellfun('isempty',regexpi(x,'granule','match'))),SHfinfo));

GCLfiles=(txtfiles(GCLidx))';
GCLfiles=cellfun(@(x) x{:}(1:end-9),GCLfiles,'UniformOutput', false);

PCLidx=logical(cellfun(@(x) sum(~cellfun('isempty',regexpi(x,'PC','match'))),SHfinfo)) | ...
    logical(cellfun(@(x) sum(~cellfun('isempty',regexpi(x,'purkinje','match'))),SHfinfo)) | ...
    logical(cellfun(@(x) sum(~cellfun('isempty',regexpi(x,'complex','match'))),SHfinfo)) | ...
    logical(cellfun(@(x) sum(~cellfun('isempty',regexpi(x,'CS','match'))),SHfinfo));

PCLfiles=(txtfiles(PCLidx))';
PCLfiles=cellfun(@(x) x{:}(1:end-9),PCLfiles,'UniformOutput', false);

MLidx=logical(cellfun(@(x) sum(~cellfun('isempty',regexpi(x,'ML','match'))),SHfinfo)) | ...
    logical(cellfun(@(x) sum(~cellfun('isempty',regexpi(x,'molecular','match'))),SHfinfo));

MLfiles=(txtfiles(MLidx))';
MLfiles=cellfun(@(x) x{:}(1:end-9),MLfiles,'UniformOutput', false);

results=[GCLfiles;PCLfiles;MLfiles];
results(1:length(GCLfiles),2)={'GCL'};
results(length(GCLfiles)+1:length(GCLfiles)+length(PCLfiles),2)={'PCL'};
results(length(GCLfiles)+length(PCLfiles)+1:end,2)={'ML'};

%%  export info
if export
    cd(directory);
    for SHfn=1:length(results)
        if strcmp(results{SHfn,1}(1),'R')
            monknum=1;
        elseif strcmp(results{SHfn,1}(1),'S')
            monknum=2;
        elseif strcmp(results{SHfn,1}(1),'H')
            monknum=3;
        end
        % get number of row in "database"
        exl = actxserver('excel.application');
        exlWkbk = exl.Workbooks;
        exlFile = exlWkbk.Open([directory 'procdata.xlsx']);
        exlSheet = exlFile.Sheets.Item(monknum);% e.g.: 2 = Sixx
        robj = exlSheet.Columns.End(4);
        numrows = robj.row;
        if numrows==1048576 %empty document
            numrows=1;
        end
        Quit(exl);
        
        [~,pfilelist] = xlsread('procdata.xlsx',monknum,['A2:A' num2str(numrows)]);
        
        if logical(sum(ismember(pfilelist,results{SHfn,1}))) %file has to have been processed already, but if for some reason not, then do not go further
            wline=find(ismember(pfilelist,results{SHfn,1}))+1;
        else
            continue
        end
        xlswrite('procdata.xlsx', {results{SHfn,2}}, monknum, sprintf('I%d',wline));
    end
end
end