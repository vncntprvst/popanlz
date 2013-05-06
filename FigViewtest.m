function varargout = FigViewtest(varargin)
% FIGVIEWTEST MATLAB code for FigViewtest.fig
%      Designed to make Rex data analysis easier, more efficient, and nicer to
%      the eye. VP 02/2012
%
%      FIGVIEWTEST, by itself, creates a new FIGVIEWTEST or raises the existing
%      singleton*.
%
%      H = FIGVIEWTEST returns the handle to a new FIGVIEWTEST or the handle to
%      the existing singleton*.
%
%      FIGVIEWTEST('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FIGVIEWTEST.M with the given input arguments.
%
%      FIGVIEWTEST('Property','Value',...) creates a new FIGVIEWTEST or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FigViewtest_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FigViewtest_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Last Modified by GUIDE v2.5 06-May-2013 11:17:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @FigViewtest_OpeningFcn, ...
    'gui_OutputFcn',  @FigViewtest_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before FigViewtest is made visible.
function FigViewtest_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FigViewtest (see VARARGIN)
global directory slash;
if strcmp(getenv('username'),'SommerVD') ||...
        strcmp(getenv('username'),'LabV') || ...
        strcmp(getenv('username'),'Purkinje')|| ...
        strcmp(getenv('username'),'vp35')
    directory = 'C:\Data\Recordings\';
elseif strcmp(getenv('username'),'DangerZone')
    directory = 'E:\data\Recordings\';
elseif strcmp(getenv('username'),'Radu')
        directory = 'E:\Spike_Sorting\';
else
    directory = 'B:\data\Recordings\';
end
slash = '\';

% Choose default command line output for FigViewtest
handles.output = hObject;

set(hObject,'DefaultTextFontName','Calibri'); %'Color',[0.9 .9 .8]

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = FigViewtest_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in arrowbackw.
function arrowbackw_Callback(hObject, eventdata, handles)
% hObject    handle to arrowbackw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% global replacespikes;
% 
% showwtrial= get(get(findobj('Tag','showtrials'),'SelectedObject'),'Tag');%selected button's tag
% if strcmp(showwtrial,'showalltrials')
%     trialnumber=str2num(get(findobj('Tag','trialnumbdisplay'),'String'))-1;
%     if trialnumber<1
%         trialnumber=str2num(get(findobj('Tag','nboftrialsdisplay'),'String'));
%     end
% elseif strcmp(showwtrial,'showoutlierstrials')
%     outliers=str2num(get(findobj('Tag','outliertrialnb'),'String'));
%     currtrial=str2num(get(findobj('Tag','trialnumbdisplay'),'String'));
%     if find(outliers==currtrial)==1
%         trialnumber=max(outliers);
%     else
%         trialnumber=outliers(find(outliers==currtrial)-1);
%     end
% elseif strcmp(showwtrial,'showbadtrials') %another day
% elseif strcmp(showwtrial,'showgoodtrials') %another day
% end
% if ~cellfun('isempty',regexp(get(findobj('Tag','filenamedisplay'),'String'),'\d\d$', 'match')) %'native' filename, without _Sp2 or _REX appendice
%     if strcmp(showwtrial,'showoutlierstrials') %was just processed
%         if replacespikes
%             filename=[get(findobj('Tag','filenamedisplay'),'String') '_Sp2'];
%         else
%             filename=[get(findobj('Tag','filenamedisplay'),'String') '_REX'];
%         end
%     else
%         filename=get(findobj('Tag','filenamedisplay'),'String')
%     end
% end
% set(findobj('Tag','trialnumbdisplay'),'String',num2str(trialnumber));
% rdd_trialdata(filename, trialnumber);

% --- Executes on button press in arrowforw.
function arrowforw_Callback(hObject, eventdata, handles)
% hObject    handle to arrowforw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% global replacespikes;
% 
% showwtrial= get(get(findobj('Tag','showtrials'),'SelectedObject'),'Tag');%selected button's tag
% if strcmp(showwtrial,'showalltrials')
%     trialnumber=str2num(get(findobj('Tag','trialnumbdisplay'),'String'))+1;
%     if trialnumber>str2num(get(findobj('Tag','nboftrialsdisplay'),'String'))
%         trialnumber=1;
%     end
% elseif strcmp(showwtrial,'showoutlierstrials')
%     outliers=str2num(get(findobj('Tag','outliertrialnb'),'String'));
%     currtrial=str2num(get(findobj('Tag','trialnumbdisplay'),'String'));
%     if find(outliers==currtrial)==length(outliers)
%         trialnumber=min(outliers);
%     else
%         trialnumber=outliers(find(outliers==currtrial)+1);
%     end
% elseif strcmp(showwtrial,'showbadtrials') %another day
% elseif strcmp(showwtrial,'showgoodtrials') %another day
% end
% if ~cellfun('isempty',regexp(get(findobj('Tag','filenamedisplay'),'String'),'\d\d$', 'match')) %'native' filename, without _Sp2 or _REX appendice
%     if strcmp(showwtrial,'showoutlierstrials') %was just processed
%         if replacespikes
%             filename=[get(findobj('Tag','filenamedisplay'),'String') '_Sp2'];
%         else
%             filename=[get(findobj('Tag','filenamedisplay'),'String') '_REX'];
%         end
%     else
%         filename=get(findobj('Tag','filenamedisplay'),'String')
%     end
% end
% set(findobj('Tag','trialnumbdisplay'),'String',num2str(trialnumber));
% rdd_trialdata(filename, trialnumber);

% --- Executes on selection change in SacDirToDisplay.
function SacDirToDisplay_Callback(hObject, eventdata, handles)
% hObject    handle to SacDirToDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SacDirToDisplay contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SacDirToDisplay


% --- Executes during object creation, after setting all properties.
function SacDirToDisplay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SacDirToDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in trueeffect.
function trueeffect_Callback(hObject, eventdata, handles)
% hObject    handle to trueeffect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global directory slash
% Matlab to Excel functions:
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/298622
% read the following threads:
% 1. <http://newsreader.mathworks.com/WebX/.ef55946>
% 2. <http://newsreader.mathworks.com/WebX/.ef549db>
% 
% The basic method is to start recording a new macro in Excel, do
% whatever you need, stop recording, then edit the macro (Alt-F8) to
% extract the relevant VB commands and matlabize-them. It's pretty
% straight forward.

%turn excel cell content to bold
h = actxserver('Excel.Application')
h.visible = true
h.Workbooks.Add
sheet = h.ActiveWorkbook.Sheets.Item(1);
sheet.Activate
cells = h.ActiveSheet.Range('C1')
set(cells.Font, 'Bold', true)

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over trueeffect.
function trueeffect_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to trueeffect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in savechg.
function savechg_Callback(hObject, eventdata, handles)
% hObject    handle to savechg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on selection (double click) of file in listbox.
% this the function that loads data,  puts name in UserData, and display
% first trial

function displayfigs_Callback(hObject, eventdata, handles)
% hObject    handle to displayfigs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global directory slash

% subjselected=get(get(findobj('Tag','monkeyselect'),'SelectedObject'),'Tag');
% taskselected=get(get(findobj('Tag','taskpanel'),'SelectedObject'),'Tag');
% evtcatselected=get(get(findobj('Tag','displayfigcat'),'SelectedObject'),'Tag');
% assumedfxonly=get(findobj('Tag','assumedeffect'),'Value');
% confirmedfxonly=get(findobj('Tag','confirmedeffect'),'Value');

% if strcmp(monkeydirselected,'rigelselect')
%     subjdir = [directory,'Rigel',slash]; %'B:\data\Recordings\Rigel';
%     procdir = [directory,'processed',slash,'Rigel',slash];
% elseif strcmp(monkeydirselected,'sixxselect')
%     subjdir = [directory,'Sixx',slash]; %'B:\data\Recordings\Sixx';
%     procdir = [directory,'processed',slash,'Sixx',slash];
% elseif strcmp(monkeydirselected,'hildaselect')
%     subjdir = [directory,'Hilda',slash]; %'B:\data\Recordings\Hilda';
%     procdir = [directory,'processed',slash,'Hilda',slash];
% elseif strcmp(monkeydirselected,'shufflesselect')
%     subjdir = [directory,'Shuffles',slash]; %'B:\data\Recordings\Hilda';
%     procdir = [directory,'processed',slash,'Shuffles',slash];
% end

if strcmp(get(gcf,'SelectionType'),'normal') && ~strcmp(eventdata,'rightclkevt') % if simple click, just higlight it, don't open
%     if get(findobj('Tag','displayfbt_files'),'Value') % works only with individual file display
%         %set uimenu content for following rightclick
%         %SelectionType 'Alternate' (Right click) doesn't work with listbox
%         %dispmenu=(get(hObject,'UIContextMenu'));
%         set(findobj('tag','displaymfiles'),'Max',1)
%         clear rightclkevt;
%         clear global tasktype;
%         listboxcontextmenu=uicontextmenu;
%         processedrexfiles = cellstr(get(hObject,'String')); % returns displaymfiles contents as cell array
%         rclk_filename = processedrexfiles{get(hObject,'Value')}; %returns selected item from displaymfiles
%         filecontent=matfile([procdir,rclk_filename,'.mat']);
%         filecodes=filecontent.allcodes;
%         curtasktype=taskdetect(filecodes);
%         if iscell(curtasktype)
%             curtasktype=cell2mat(curtasktype);
%         end
%         disptask=uimenu('Parent',listboxcontextmenu,'Label',curtasktype);
%         set(hObject,'UIContextMenu',listboxcontextmenu);
%     else
%         set(hObject,'Max',2);
%         rightclkevt='rightclkevt';
%         multianalyzepopup=uicontextmenu();
%         set(hObject,'uicontextmenu',multianalyzepopup);
%         item1 = uimenu(multianalyzepopup,'Label','Analyze','Callback',{@(src,evt,handles) displaymfiles_Callback(hObject,rightclkevt,handles),hObject});  %disp(get(hObject,'Value'))
%         %displaymfiles_Callback(hObject, eventdata, handles)
%     end
elseif strcmp(get(gcf,'SelectionType'),'open') %|| strcmp(eventdata,'rightclkevt')
    cd(directory);

    dirfignames = cellstr(get(hObject,'String')); % returns displaymfiles contents as cell array
    selectionnm = {dirfignames{get(hObject,'Value')}}; %returns selected item from displaymfiles
    openimageinaxes(selectionnm,handles.figview);
    
    % write file name in the file name display
    set(findobj('Tag','filenamedisplay'),'String',regexpi(selectionnm{:},'[a-z]\d+[a-z]\d[a-z]\d_\d+','match'));
    % same thing for task
    rectasklims=strfind(selectionnm{:},'_');
    set(findobj('Tag','taskdisplay'),'String',selectionnm{:}(rectasklims(3)+1:rectasklims(4)-1));    
end


% --- Executes during object creation, after setting all properties.
function displayfigs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to displayfigs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global directory slash

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% determines computer type
archst  = computer('arch');

if strcmp(archst, 'maci64')
    name = getenv('USER');
    if strcmp(name, 'nick')
        directory = '/Users/nick/Dropbox/filesforNick/';
    elseif strcmp(name, 'Frank')
        directory = '/Users/Frank/Desktop/monkeylab/data/';
    elseif strcmp(name, 'zacharyabzug')
        directory = '/Users/zacharyabzug/Desktop/zackdata/';
    end
    slash = '/';
elseif strcmp(archst, 'win32') || strcmp(archst, 'win64')
    if strcmp(getenv('username'),'SommerVD') || ...
            strcmp(getenv('username'),'LabV') || ...
            strcmp(getenv('username'),'Purkinje') || ...
            strcmp(getenv('username'),'vp35')
        directory = 'C:\Data\Recordings\';
    elseif strcmp(getenv('username'),'DangerZone')
        directory = 'E:\data\Recordings\';
    elseif strcmp(getenv('username'),'Radu')
        directory = 'E:\Spike_Sorting\';
    else
        directory = 'B:\data\Recordings\';
    end
    slash = '\';
end

%listing figures directory
figdir{1} = dir([directory,'figures',slash,'sac',slash]);
figdir{2} = dir([directory,'figures',slash,'vis',slash]);
figdir{3} = dir([directory,'figures',slash,'rew',slash]);


listcat=get(get(findobj('Tag','displayfigcat'),'SelectedObject'),'Tag');
if strcmp(listcat,'displaysacfigs')
    figcat=1;
elseif strcmp(listcat,'displayvisfigs')
    figcat=2;
elseif strcmp(listcat,'displayrewfigs')
    figcat=3;
elseif strcmp(listcat, 'displayallfigs')
    figcat=0;
end

%%  list figures
if figcat>0
   dirfignames={figdir{figcat}.name};
   % Order by date
   filedates=[figdir{figcat}.datenum];
else
    dirfignames={figdir{:}.name}; %incorrect
    filedates=cell2mat({figdir(:).datenum});
end

[filedates,fdateidx] = sort(filedates,'descend');
dirfignames = dirfignames(fdateidx);
dirfignames = cellfun(@(x) x(1:end-4), dirfignames, 'UniformOutput', false);
dirfignames=dirfignames(~cellfun('isempty',regexpi(dirfignames,'^\w','match')));

%only show files for selected subject
subjdir= get(get(findobj('Tag','monkeyselect'),'SelectedObject'),'Tag');
if strcmp(subjdir,'rigelselect')
    monknum=1;
    subjinitial='H';
elseif strcmp(subjdir,'sixxselect')
    monknum=2;
    subjinitial='S';
elseif strcmp(subjdir,'hildaselect')
    monknum=3;
    subjinitial='H';
elseif strcmp(subjdir, 'shufflesselect')
    monknum=4;
    subjinitial='Sh';
elseif strcmp(subjdir, 'allsubjects')
    monknum=0;
end

if monknum>0
    dirfignames=dirfignames(~cellfun('isempty',regexpi(dirfignames,['^',subjinitial],'match')));
end

%only show files for selected task
taskselected=get(get(findobj('Tag','taskpanel'),'SelectedObject'),'Tag');
if strcmp(taskselected,'stsacselect')
    figtask='st_saccades';
elseif strcmp(taskselected,'cmdselect')
    figtask='gapstop';
elseif strcmp(taskselected, 'tokselect')
    figtask='tokens';
elseif strcmp(taskselected, 'optilocselect')
    figtask='optiloc';
elseif strcmp(taskselected,'alltaskselect')
    figtask='all';
end

if ~strcmp(figtask,'all')
    dirfignames=dirfignames(~cellfun('isempty',regexpi(dirfignames,figtask,'match')));
end

%% get activity data from process file

            cd(directory);
            % get current processed file list and activity from excel file
            for monknum=1:3
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
            
                [~,pfilelist{monknum}] = xlsread('procdata.xlsx',monknum,['A2:A' num2str(numrows)]);
                [activlist{monknum}] = xlsread('procdata.xlsx',monknum,['K2:K' num2str(numrows)]);
                if length(activlist{monknum})<length(pfilelist{monknum})
                    activlist{monknum}=[activlist{monknum};nan(length(pfilelist{monknum})-length(activlist{monknum}),1)];
                end
                [~,activetype{monknum}] = xlsread('procdata.xlsx',monknum,['L2:L' num2str(numrows)]);
                [~,roughloc{monknum}] = xlsread('procdata.xlsx',monknum,['H2:H' num2str(numrows)]);
            end
            % get FontWeight from excel sheet. Not sure how to do. Check
            % this page for ideas
%             http://www.mathworks.com/matlabcentral/fileexchange/4981-xlsfont-xlsalign-xlsborder-xlswordart-xlscomment/content/xlsfont.m
%           maybe try: get(Excel.Selection.Font,'FontStyle');

if get(findobj('Tag','assumedeffect'),'Value'); % keep only files with assumed effect
origfilename=regexpi(dirfignames,'[a-z]\d+[a-z]\d[a-z]\d_\d+','match');
fullpfilelist=cat(1,pfilelist{:});
fullactivelist=cat(1,activlist{:});
activefiles=fullpfilelist(fullactivelist>1);
dirfignames=dirfignames(ismember([origfilename{:}],activefiles));
% corresflist=cellfun(@(x) ismember([origfilename{:}],x),pfilelist,'UniformOutput',false);
% dirfignames=dirfignames(fullactivelist(cat(1,corresflist{:}))>1);
end

%confirmedfxonly=get(findobj('Tag','confirmedeffect'),'Value');

% create structure to pass the activity data
initdata = {struct('pfilelist',pfilelist,'activlist',activlist,'activetype',activetype,'roughloc',roughloc)};

set(hObject,'string',dirfignames);
guidata(findobj('tag','Effectselect'),initdata);
guidata(hObject,figdir);


% --- Executes when selected object is changed in displayfigcat.
function displayfigcat_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in displayfigcat
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global directory slash;

% set(eventdata.OldValue, 'BackgroundColor', [0.9608    0.9216    0.9216]);

if get(findobj('Tag','rigelselect'),'Value')
    dirlisting = dir([directory,'processed',slash,'Rigel',slash]); %('B:\data\Recordings\processed\Rigel');
elseif get(findobj('Tag','sixxselect'),'Value')
    dirlisting = dir([directory,'processed',slash,'Sixx',slash]); %('B:\data\Recordings\processed\Sixx');\
elseif get(findobj('Tag','hildaselect'),'Value')
    dirlisting = dir([directory,'processed',slash,'Hilda',slash]); %('B:\data\Recordings\processed\Sixx');
elseif get(findobj('Tag','shufflesselect'),'Value')
    dirlisting = dir([directory,'processed',slash,'Shuffles',slash]); %('B:\data\Recordings\processed\Sixx');
end
fileNames = {dirlisting.name};  % Put the file names in a cell array

if hObject==findobj('Tag','displayfbt_files')
    
    % Order by date
    filedates=cell2mat({dirlisting(:).datenum});
    [filedates,fdateidx] = sort(filedates,'descend');
    dirlisting = {dirlisting(:).name};
    dirlisting = dirlisting(fdateidx);
    dirlisting = dirlisting(~cellfun('isempty',strfind(dirlisting,'mat')));
    dirlisting = dirlisting(cellfun('isempty',strfind(dirlisting,'myBreakpoints')));
    dirlisting = cellfun(@(x) x(1:end-4), dirlisting, 'UniformOutput', false);
    set(findobj('Tag','displaymfiles'),'string',dirlisting);
    
elseif hObject==findobj('Tag','displayfbt_session')
    
    if get(findobj('Tag','rigelselect'),'Value')
        index = regexpi(fileNames,...              % Match a file name if it begins
            '^R\d+','match');           % with the letter 'R' followed by a set of digits 1 or larger
        inFiles = index(~cellfun(@isempty,index));  % Get the names of the matching files in a cell array
        sessionNumbers = cellfun(@(x) strrep(x, 'R', ' '), inFiles, 'UniformOutput', false);
    elseif get(findobj('Tag','sixxselect'),'Value')
        index = regexpi(fileNames,...              % Match a file name if it begins
            '^S\d+', 'match');           % with the letter 'S' followed by a set of digits 1 or larger
        inFiles = index(~cellfun(@isempty,index));  % Get the names of the matching files in a cell array
        sessionNumbers = cellfun(@(x) strrep(x, 'S', ' '), inFiles, 'UniformOutput', false);
    elseif get(findobj('Tag','hildaselect'),'Value')
        index = regexpi(fileNames,...              % Match a file name if it begins
            '^H\d+', 'match');           % with the letter 'S' followed by a set of digits 1 or larger
        inFiles = index(~cellfun(@isempty,index));  % Get the names of the matching files in a cell array
        sessionNumbers = cellfun(@(x) strrep(x, 'H', ' '), inFiles, 'UniformOutput', false);
    elseif get(findobj('Tag','shufflesselect'),'Value')
        index = regexpi(fileNames,...              % Match a file name if it begins
            '^S\d+', 'match');           % with the letter 'S' followed by a set of digits 1 or larger
        inFiles = index(~cellfun(@isempty,index));  % Get the names of the matching files in a cell array
        sessionNumbers = cellfun(@(x) strrep(x, 'S', ' '), inFiles, 'UniformOutput', false);
    end
    
    if ~isempty(sessionNumbers)
        dispsession = cat(1,sessionNumbers{:});
        dispsession = unique(dispsession); % finds unique session numbers
        [~,sessionidx]=sort(str2double(dispsession),'descend');
        dispsession = dispsession(sessionidx); %sort cell arrays descending
        dispsession = strcat('Session', dispsession);
        set(findobj('Tag','displaymfiles'),'string', dispsession);
    else
        set(findobj('Tag','displaymfiles'),'string','');
    end
    
elseif hObject==findobj('Tag','displayfbt_grid')
    
    index = regexpi(fileNames,...              % Match a file name if it begins
        '[a-z]\d[a-z]\d','match');           % with the letter 'R' followed by a set of digits 1 or larger
    gridLocations = index(~cellfun(@isempty,index));  % Get the names of the matching files in a cell array
    
    if ~isempty(gridLocations)
        displocation = cat(1,gridLocations{:});
        displocation = unique(displocation); % finds unique session numbers
        [~,sessionidx]=sort(str2double(displocation),'descend');
        displocation = displocation(sessionidx); %sort cell arrays descending
        set(findobj('Tag','displaymfiles'),'string', displocation);
    else
        set(findobj('Tag','displaymfiles'),'string','');
    end
    
end

% --- Executes on button press in unconfirmedeffect.
function unconfirmedeffect_Callback(hObject, eventdata, handles)
% hObject    handle to unconfirmedeffect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in confirmedeffect.
function confirmedeffect_Callback(hObject, eventdata, handles)
% hObject    handle to confirmedeffect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of confirmedeffect


% --- Executes on button press in assumedeffect.
function assumedeffect_Callback(hObject, eventdata, handles)
% hObject    handle to assumedeffect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of assumedeffect


% --- Executes on key press with focus on displayvisfigs and none of its controls.
function displayvisfigs_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to displayvisfigs (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function sixxselect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sixxselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function hildaselect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hildaselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function taskpanel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to taskpanel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function filenamedisplay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filenamedisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function monkeyselect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to monkeyselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function taskdisplay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to taskdisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
