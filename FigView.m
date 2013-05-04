function varargout = FigView(varargin)
% FIGVIEW M-file for FigView.fig
%      FIGVIEW, by itself, creates a new FIGVIEW or raises the existing
%      singleton*.
%
%      H = FIGVIEW returns the handle to a new FIGVIEW or the handle to
%      the existing singleton*.
%
%      FIGVIEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FIGVIEW.M with the given input arguments.
%
%      FIGVIEW('Property','Value',...) creates a new FIGVIEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FigView_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FigView_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FigView

% Last Modified by GUIDE v2.5 03-May-2013 19:53:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FigView_OpeningFcn, ...
                   'gui_OutputFcn',  @FigView_OutputFcn, ...
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


% --- Executes just before FigView is made visible.
function FigView_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FigView (see VARARGIN)

% Choose default command line output for FigView
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FigView wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FigView_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
