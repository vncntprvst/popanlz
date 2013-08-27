
if strcmp(getenv('username'),'SommerVD') || strcmp(getenv('username'),'vp35')
    directory = 'C:\Data\Recordings\';
elseif  strcmp(getenv('username'),'DangerZone')
    directory = 'E:\Data\Recordings\';
else
    directory = 'B:\data\Recordings\';
end
slash = '\';

%% get gapstop filenames from cmd excel file


%Get number rows
exl = actxserver('excel.application');
exlWkbk = exl.Workbooks;
exlFile = exlWkbk.Open([directory 'cmd.xlsx']);
exlSheet = exlFile.Sheets.Item(1);% all selected data in first sheet
robj = exlSheet.Columns.End(4);
numrows = robj.row;
% if numrows==1048576 %empty document
%     numrows=1;
% end
Quit(exl);

% read A (File Name)
[~,filestoprint] = xlsread([directory 'cmd.xlsx'],1,['A2:A' num2str(numrows)]);

cmdfigdir = [directory,'figures',slash,'cmd',slash];
cmdfigdirlisting=dir(cmdfigdir);
cmdfigdirfileNames={cmdfigdirlisting.name};

for filenum=1:length(filestoprint)
    
    if logical(sum(~cellfun('isempty',regexpi(cmdfigdirfileNames,[filestoprint{filenum} '_REX_NSSvsCSS_tgt.png'],'match')))) &&...
            logical(sum(~cellfun('isempty',regexpi(cmdfigdirfileNames,[filestoprint{filenum} '_REX_NSSvsNCSS_sac.png'],'match'))))
        
        continue
    
    else
        if strcmp('R',filestoprint{filenum}(1))
            subject='Rigel';
            procdir = [directory,'processed',slash,'Rigel',slash];
        elseif strcmp('S',filestoprint{filenum}(1))
            subject='Sixx';
            procdir = [directory,'processed',slash,'Sixx',slash];
        elseif strcmp('H',filestoprint{filenum}(1))
            subject='Hilda';
            procdir = [directory,'processed',slash,'Hilda',slash];
        end
        
        procdirlisting=dir(procdir);
        procdirfileNames={procdirlisting.name};
        procdirfileNames = regexprep(procdirfileNames, '(_REX.mat$)|(_Sp2.mat$)','');
        loadfile=procdirfileNames(~cellfun('isempty',regexpi(procdirfileNames,filestoprint{filenum},'match')));
        % if two versions (REX and Spike2), chose Spike2
        if length(loadfile)>1
            if ~isempty(loadfile(~cellfun('isempty',regexpi(loadfile,'Sp2','match'))))
                loadfile=loadfile(~cellfun('isempty',regexpi(loadfile,'Sp2','match')));
            else
                loadfile={[filestoprint{filenum} '_REX']};
            end   
        end
        
        [~, trialdirs] = data_info(loadfile{:}, 1, 1); %reload file: yes (shouldn't happen, though), skip unprocessed files: yes
        getaligndata{filenum} = rdd_rasters_sdf(loadfile{:}, trialdirs, 0); % align data, don't plot rasters
    end  
end



