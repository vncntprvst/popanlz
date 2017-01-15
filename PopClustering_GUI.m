function varargout = PopClustering_GUI(varargin)
% POPCLUSTERING_GUI MATLAB code for PopClustering_GUI.fig
%      POPCLUSTERING_GUI, by itself, creates a new POPCLUSTERING_GUI or raises the existing
%      singleton*.
%
%      H = POPCLUSTERING_GUI returns the handle to a new POPCLUSTERING_GUI or the handle to
%      the existing singleton*.
%
%      POPCLUSTERING_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in POPCLUSTERING_GUI.M with the given input arguments.
%
%      POPCLUSTERING_GUI('Property','Value',...) creates a new POPCLUSTERING_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PopClustering_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PopClustering_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%       
%       Examples:   PopClustering_GUI % no data set selected
%                   PopClustering_GUI('cDn_gsdata') %loads cDn_gsdata
%                   PopClustering_GUI(data,dataset) %uses input clustering data
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PopClustering_GUI

% Last Modified by GUIDE v2.5 30-Dec-2015 14:43:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PopClustering_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @PopClustering_GUI_OutputFcn, ...
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


% --- Executes just before PopClustering_GUI is made visible.
function PopClustering_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% varargin   command line arguments to PopClustering_GUI (see VARARGIN)

%% settings
handles.userinfo=SetUserDir;

cd(handles.userinfo.syncdir);
if isempty(varargin) || size(varargin,2)==1
    if ~isempty(varargin)
        handles.dataset=varargin{1};
    else
        % select dataset to load
        datasets={'cDn_gsdata','top_cortex_gsdata'};
        [datasetNum,okgo] = listdlg('PromptString','Select dataset to load:',...
            'SelectionMode','single',...
            'ListString',datasets,...
            'OKString','Load');
        if okgo
            handles.dataset=datasets{datasetNum};
        else %default
            handles.dataset='cDn_gsdata';
        end
    end
    data=load([handles.dataset '.mat']); %cDn_gsdata.mat  top_cortex_gsdata.mat
    dataField=cell2mat(fieldnames(data));
    handles.unitsDBinfo=data.(dataField).alldb;
    handles.dbConn = connect2DB('vp_sldata');
    % get cluster index from database
    unitList=cellfun(@(x) x.unit_id, handles.unitsDBinfo);
    handles.clusterIdx=fetch(handles.dbConn,['SELECT profile_type FROM clusters c WHERE c.cluster_id IN (' ...
        sprintf('%.0f,' ,unitList(1:end-1)) num2str(unitList(end)) ')']);
    [~,resort]=sort(unitList);[~,resort]=sort(resort);handles.clusterIdx=[handles.clusterIdx{resort}]';
    %remove bad apples
    handles.unitsDBinfo=handles.unitsDBinfo(~isnan(handles.clusterIdx));
    handles.clusterIdx=handles.clusterIdx(~isnan(handles.clusterIdx));
    % process data
%     [handles.unitsProfile.sacresps,handles.unitsProfile.bnorm_sacresps,handles.unitsProfile.rnorm_sacresps]=comp_sacresp(data,dataField);
    % normalize data
    data.(dataField).normData(:,1),);% send sac epoch data only (cut each raster? )
    handles.unitsProfile.norm_sacresps=RespNormalization(data, data.(dataField).normFactor)
    

else
    handles.dataset=varargin{1};
    data=varargin{2};
    %populate handles
    handles.unitsProfile=data.cellprofile;
    handles.clusterIdx=data.clusteridx;
	handles.unitsDBinfo=data.dbinfo;
    handles.dbConn=data.dbconn;
end

% get list of profiles
unitList=cellfun(@(x) x.unit_id, handles.unitsDBinfo);
handles.unitsProfileType=fetch(handles.dbConn,['SELECT profile FROM clusters c WHERE c.cluster_id IN (' ...
        sprintf('%.0f,' ,unitList(1:end-1)) num2str(unitList(end)) ')']);
[~,resort]=sort(unitList);[~,resort]=sort(resort);handles.unitsProfileType=handles.unitsProfileType(resort);

% Choose default command line output for PopClustering_GUI
%% initial plots 
handles.output = hObject;
handles.fcm=colormap(hObject,lines);
% plot average waveforms
% cmrtitles={'Unsorted','early peak','ramp & fall','burst','ramp all the way'};
clusid=unique(handles.clusterIdx);
clusterProfiles=cell(length(unique(handles.clusterIdx)),1);
for mclussp=1:length(unique(handles.clusterIdx))
%      %findobj('Tag',['cluster' num2str(mclussp)])
    clusterTag=['cluster' num2str(mclussp)];
%     % add waveform to handle structure
    handles.([clusterTag '_wf'])= nanmean(handles.unitsProfile.bnorm_sacresps(handles.clusterIdx==clusid(mclussp),:));
    % collect profile names
    [DBprofiles,~,DBprofilesIdx]=unique(handles.unitsProfileType(handles.clusterIdx==clusid(mclussp),:));
    handles.([clusterTag '_Profile'])= DBprofiles{mode(DBprofilesIdx)};
    clusterProfiles{mclussp}=DBprofiles{mode(DBprofilesIdx)};
    axes(handles.(clusterTag));
%     handles.(clusterTag)=
    plot(handles.(clusterTag),handles.([clusterTag '_wf']),'Color',handles.fcm(mclussp,:),'LineWidth',2);
%      text(100,mean(get(gca,'ylim')),['n=' num2str(sum(handles.clusterIdx==clusid(mclussp)))]);
%     xlabel('Time')
%     ylabel('Norm. Firing rate')
% %     title(cmrtitles{mclussp})
    title([clusterTag ' (' num2str(clusid(mclussp)) '), n=' num2str(sum(handles.clusterIdx==clusid(mclussp)))])
    set(gca,'xtick',1:100:max(get(gca,'xlim')),'xticklabel',-(max(get(gca,'xlim'))/2):100:max(get(gca,'xlim'))/2,'TickDir','out');
    set(gca,'Color','white','FontSize',8,'FontName','calibri');
    axis(gca,'tight'); box off;
end
% set(hObject,'visible','off')

%populate cluster number listbox
handles.currentUnitDisplay=1;
set(handles.changeclustype_listbox,'String',num2str(clusid));
set(handles.changeclustype_listbox,'value',find(clusid==handles.clusterIdx(handles.currentUnitDisplay)));

%populate cluster profiles names box
set(handles.NameClusters_box,'String',clusterProfiles);

%plot first unit and display info
axes(handles.mainplot);
plot(handles.mainplot,handles.unitsProfile.bnorm_sacresps(handles.currentUnitDisplay,:),'Color',...
    handles.fcm(clusid==handles.clusterIdx(handles.currentUnitDisplay),:),'LineWidth',2);
    set(gca,'xtick',1:100:max(get(gca,'xlim')),'xticklabel',-(max(get(gca,'xlim'))/2):100:max(get(gca,'xlim'))/2,'TickDir','out');
    set(gca,'Color','white','FontSize',8,'FontName','calibri');
    axis(gca,'tight'); box off;
recName=fetch(handles.dbConn,['SELECT r.e_file FROM recordings r WHERE r.recording_id = ' num2str(handles.unitsDBinfo{handles.currentUnitDisplay, 1}.rec_id) ';']);
set(handles.recname_box,'string',{recName{:}(1:end-1)},'FontSize',10,'FontName','calibri','FontWeight','bold');
set(handles.clusternumber_box,'string',{num2str(handles.clusterIdx(handles.currentUnitDisplay))},'FontSize',10,'FontName','calibri','FontWeight','bold');
set(handles.clustertype_box,'string',{handles.unitsProfileType{handles.currentUnitDisplay}},'FontSize',10,'FontName','calibri','FontWeight','bold');
set(hObject,'visible','off')
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PopClustering_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PopClustering_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in previousplot.
function previousplot_Callback(hObject, eventdata, handles)
% hObject    handle to previousplot (see GCBO)

switch handles.currentUnitDisplay
    case 1
        handles.currentUnitDisplay=size(handles.unitsDBinfo,1);
    otherwise
        handles.currentUnitDisplay=handles.currentUnitDisplay-1;
end
%update unit number display
set(handles.unitNumber,'String',num2str(handles.currentUnitDisplay));
%update plot
clusid=unique(handles.clusterIdx);
axes(handles.mainplot);
plot(handles.mainplot,handles.unitsProfile.bnorm_sacresps(handles.currentUnitDisplay,:),'Color',...
    handles.fcm(clusid==handles.clusterIdx(handles.currentUnitDisplay),:),'LineWidth',2);
    set(gca,'xtick',1:100:max(get(gca,'xlim')),'xticklabel',-(max(get(gca,'xlim'))/2):100:max(get(gca,'xlim'))/2,'TickDir','out');
    set(gca,'Color','white','FontSize',8,'FontName','calibri');
    axis(gca,'tight'); box off;
recName=fetch(handles.dbConn,['SELECT r.e_file FROM recordings r WHERE r.recording_id = ' num2str(handles.unitsDBinfo{handles.currentUnitDisplay, 1}.rec_id) ';']);
set(handles.recname_box,'string',{recName{:}(1:end-1)},'FontSize',10,'FontName','calibri','FontWeight','bold');
set(handles.clusternumber_box,'string',{num2str(handles.clusterIdx(handles.currentUnitDisplay))},'FontSize',10,'FontName','calibri','FontWeight','bold');
set(handles.clustertype_box,'string',{handles.unitsProfileType{handles.currentUnitDisplay}},'FontSize',10,'FontName','calibri','FontWeight','bold');
set(handles.changeclustype_listbox,'value',find(clusid==handles.clusterIdx(handles.currentUnitDisplay)));
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in nextplot.
function nextplot_Callback(hObject, eventdata, handles)
% hObject    handle to nextplot (see GCBO)
switch handles.currentUnitDisplay
    case size(handles.unitsDBinfo,1)
        handles.currentUnitDisplay=1;
    otherwise
        handles.currentUnitDisplay=handles.currentUnitDisplay+1;
end
%update unit number display
set(handles.unitNumber,'String',num2str(handles.currentUnitDisplay));
%update plot
clusid=unique(handles.clusterIdx);
axes(handles.mainplot);
plot(handles.mainplot,handles.unitsProfile.bnorm_sacresps(handles.currentUnitDisplay,:),'Color',...
    handles.fcm(clusid==handles.clusterIdx(handles.currentUnitDisplay),:),'LineWidth',2);
    set(gca,'xtick',1:100:max(get(gca,'xlim')),'xticklabel',-(max(get(gca,'xlim'))/2):100:max(get(gca,'xlim'))/2,'TickDir','out');
    set(gca,'Color','white','FontSize',8,'FontName','calibri');
    axis(gca,'tight'); box off;
recName=fetch(handles.dbConn,['SELECT r.e_file FROM recordings r WHERE r.recording_id = ' num2str(handles.unitsDBinfo{handles.currentUnitDisplay, 1}.rec_id) ';']);
set(handles.recname_box,'string',{recName{:}(1:end-1)},'FontSize',10,'FontName','calibri','FontWeight','bold');
set(handles.clusternumber_box,'string',{num2str(handles.clusterIdx(handles.currentUnitDisplay))},'FontSize',10,'FontName','calibri','FontWeight','bold');
set(handles.clustertype_box,'string',{handles.unitsProfileType{handles.currentUnitDisplay}},'FontSize',10,'FontName','calibri','FontWeight','bold');
set(handles.changeclustype_listbox,'value',find(clusid==handles.clusterIdx(handles.currentUnitDisplay)));
% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in changeclustype_listbox.
function changeclustype_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to changeclustype_listbox (see GCBO)

% --- Executes during object creation, after setting all properties.
function changeclustype_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to changeclustype_listbox (see GCBO)

% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in apply_changeclustype.
function apply_changeclustype_Callback(hObject, eventdata, handles)
% hObject    handle to apply_changeclustype (see GCBO)

%get value
clusterList = cellstr(get(handles.changeclustype_listbox,'String'));
originalUnitClusIdx=handles.clusterIdx(handles.currentUnitDisplay);
handles.clusterIdx(handles.currentUnitDisplay)=str2double(clusterList{get(handles.changeclustype_listbox,'Value')});
%update plot color and unit cluster info
clusid=unique(handles.clusterIdx);
axes(handles.mainplot);
plot(handles.mainplot,handles.unitsProfile.bnorm_sacresps(handles.currentUnitDisplay,:),'Color',...
    handles.fcm(clusid==handles.clusterIdx(handles.currentUnitDisplay),:),'LineWidth',2);
set(gca,'xtick',1:100:max(get(gca,'xlim')),'xticklabel',-(max(get(gca,'xlim'))/2):100:max(get(gca,'xlim'))/2,'TickDir','out');
set(gca,'Color','white','FontSize',8,'FontName','calibri');
axis(gca,'tight'); box off;
set(handles.clusternumber_box,'string',{num2str(handles.clusterIdx(handles.currentUnitDisplay))},'FontSize',10,'FontName','calibri','FontWeight','bold');
set(handles.changeclustype_listbox,'value',find(clusid==handles.clusterIdx(handles.currentUnitDisplay)));

%% update cluster profiles and 'n' values
%old one
clusNum=find(clusid==originalUnitClusIdx);
clusterTag=['cluster' num2str(clusNum)];
%change waveform in handle structure
handles.([clusterTag '_wf'])= nanmean(handles.unitsProfile.bnorm_sacresps(handles.clusterIdx==clusid(clusNum),:));
% update cluster plot
axes(handles.(clusterTag));
plot(handles.(clusterTag),handles.([clusterTag '_wf']),'Color',handles.fcm(clusNum,:),'LineWidth',2);
title([clusterTag ' (' num2str(clusid(clusNum)) '), n=' num2str(sum(handles.clusterIdx==clusid(clusNum)))])
set(gca,'xtick',1:100:max(get(gca,'xlim')),'xticklabel',-(max(get(gca,'xlim'))/2):100:max(get(gca,'xlim'))/2,'TickDir','out');
set(gca,'Color','white','FontSize',8,'FontName','calibri');
axis(gca,'tight'); box off;
% new one
clusNum=find(clusid==handles.clusterIdx(handles.currentUnitDisplay));
clusterTag=['cluster' num2str(clusNum)];
%change waveform in handle structure
handles.([clusterTag '_wf'])= nanmean(handles.unitsProfile.bnorm_sacresps(handles.clusterIdx==clusid(clusNum),:));
% update cluster plot
axes(handles.(clusterTag));
plot(handles.(clusterTag),handles.([clusterTag '_wf']),'Color',handles.fcm(clusNum,:),'LineWidth',2);
title([clusterTag ' (' num2str(clusid(clusNum)) '), n=' num2str(sum(handles.clusterIdx==clusid(clusNum)))])
set(gca,'xtick',1:100:max(get(gca,'xlim')),'xticklabel',-(max(get(gca,'xlim'))/2):100:max(get(gca,'xlim'))/2,'TickDir','out');
set(gca,'Color','white','FontSize',8,'FontName','calibri');
axis(gca,'tight'); box off;

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in commit_changestodb.
function commit_changestodb_Callback(hObject, eventdata, handles)
% hObject    handle to commit_changestodb (see GCBO)

unitList=cellfun(@(x) x.unit_id, handles.unitsDBinfo);
clustypes=get(handles.NameClusters_box,'String');
clusid=unique(handles.clusterIdx);
% update handle
for clus=1:size(clustypes,1)
handles.unitsProfileType(handles.clusterIdx==clusid(clus))=clustypes(clus);
end
% and commit
addProfile(handles.unitsProfileType, handles.clusterIdx, unitList, handles.dbConn);
% update info
set(handles.clustertype_box,'string',handles.unitsProfileType(handles.currentUnitDisplay),'FontSize',10,'FontName','calibri','FontWeight','bold');

function NameClusters_box_Callback(hObject, eventdata, handles)
% hObject    handle to NameClusters_box (see GCBO)

% Hints: get(hObject,'String') returns contents of NameClusters_box as text
%        str2double(get(hObject,'String')) returns contents of NameClusters_box as a double
% set(handles.NameClusters_box,'String',clusterProfiles);

% --- Executes during object creation, after setting all properties.
function NameClusters_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NameClusters_box (see GCBO)
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in printHeatmap.
function printHeatmap_Callback(hObject, eventdata, handles)
% hObject    handle to printHeatmap (see GCBO)

[sortedClusterIdx,sortidx]=sort(handles.clusterIdx);
clusid=unique(handles.clusterIdx);
%population raster heat map
figure('name','population raster heatmap','position',[90 450 560 570])
% subplot(1,20,1:9)
% imagesc(1:size(handles.unitsProfile.rnorm_sacresps,2),1:size(handles.unitsProfile.rnorm_sacresps,1),handles.unitsProfile.rnorm_sacresps)
% set(gca,'FontSize',15);
% xlabel('Time')
% ylabel('Neuron #')
% title('Unsorted')
heatmaph=subplot(1,10,1:9);
imagesc(1:size(handles.unitsProfile.rnorm_sacresps,2),1:size(handles.unitsProfile.rnorm_sacresps,1),handles.unitsProfile.rnorm_sacresps(sortidx(sortedClusterIdx~=clusid(1)),:))
colormap(heatmaph,parula);
set(gca,'xtick',1:100:max(floor(get(gca,'xlim'))),'xticklabel',-(max(floor(get(gca,'xlim')))/2):100:max(floor(get(gca,'xlim')))/2,'TickDir','out');
set(gca,'Color','white','FontSize',14,'FontName','calibri');
axis(gca,'tight'); box off;
xlabel('Time')
ylabel('Neuron #')
title('Unit Response - normalized')
clusbarh=subplot(1,10,10);
clusrange=false(size(sortedClusterIdx(sortedClusterIdx~=clusid(1)),1),1);
for clusNum=2:2:length(clusid)
    clusrange=clusrange | sortedClusterIdx(sortedClusterIdx~=clusid(1))==clusid(clusNum);
end
clusrange=repmat(single(clusrange),1,5);
imagesc(clusrange);
colormap(clusbarh,bone)
set(gca,'XTick', [],'XTickLabel',[],'YTick', [],'YTickLabel',[]);

% export figure
exportfigname=[handles.userinfo.syncdir(1:regexp(handles.userinfo.syncdir,'\w+$')-1) 'figures\' handles.dataset '_heatmap'];
%print png
print(gcf, '-dpng', '-noui', '-opengl','-r600', exportfigname);
%print svg
print(gcf, '-dsvg', '-noui', exportfigname);
% plot2svg([exportfigname,'.svg'],gcf, 'png');
% close(gcf)

% --- Executes on button press in printProfiles.
function printProfiles_Callback(hObject, eventdata, handles)
% hObject    handle to printProfiles (see GCBO)

clustypes=get(handles.NameClusters_box,'String');
clusid=unique(handles.clusterIdx);
% plot mean clusters profile
figure('name','clusters mean response','position',[1120 -100 540 760])
for clusNum=2:length(clusid) 
    subplot(length(clusid)-1,1,clusNum-1)
    clusterTag=['cluster' num2str(clusNum)];
    plot(zscore(handles.([clusterTag '_wf']),0,2),'Color',handles.fcm(clusNum,:),'LineWidth',2);
    title([clustypes{clusNum} ', n=' num2str(sum(handles.clusterIdx==clusid(clusNum)))])
    set(gca,'xtick',1:100:max(get(gca,'xlim')),'xticklabel',-(max(get(gca,'xlim'))/2):100:max(get(gca,'xlim'))/2,'TickDir','out');
    set(gca,'ylim',[-2 3]);
    set(gca,'Color','white','FontSize',8,'FontName','calibri');
%     axis(gca,'tight'); 
    box off;
end
% export figure
exportfigname=[handles.userinfo.syncdir(1:regexp(handles.userinfo.syncdir,'\w+$')-1)...
    'figures\' handles.dataset '_cluster_profile'];
%print png
print(gcf, '-dpng', '-noui', '-opengl','-r600', exportfigname);
%print svg
print(gcf, '-dsvg', '-noui', exportfigname);
% plot2svg([exportfigname,'.svg'],gcf, 'png');
% close(gcf)