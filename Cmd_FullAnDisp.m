%% directory
global directory slash;
if strcmp(getenv('username'),'SommerVD') || strcmp(getenv('username'),'vp35')
    directory = 'C:\Data\Recordings\';
elseif  strcmp(getenv('username'),'DangerZone')
    directory = 'E:\Data\Recordings\';
else
    directory = 'B:\data\Recordings\';
end
slash = '\';
%% files
%dentate
CmdFileName={'S113L4A5_13500';'S114L4A5_14321';'S122L4A5_14010';'S126A4L6_13690';...
    'R132L4P4_20152';'S121L4A5_13091';'S112l4a5_12971';...
    'S117L4A6_12741';'S118L4A5_13081';'S115L4A6_12871';'S116L4A6_15431';...
    'H62L5A5_22211';'S96L2A5_13651';'S99L2A5_13242';'H56L5A5_21502';...
    'S116L4A6_15450';'H25L5A5_22760';'H53L5A5_20901';'S123L4A6_13691';...
    'S125L4A6_13990';'R167L1A3_20671';'R166L1A3_20371';'S120L4A5_12912';...
    'S114L4A5_13650';'S111L4A5_14392';'S119L4A5_14391';'S105L1A8_16050';...
    'S101L2A5_11251';'S89L4A5_12401';'S102L3A4_13001';'S111L4A5_13502';...
    'H56L5A5_21321';'R97L7A1_19001'};
% 
% 
% CmdFileName={''};
%For cancellation, for the moment keep 
% CmdFileName={'S96L2A5_13651';'S115L4A6_12871'};
%+ try to fix 'H53L5A5_20901'

%For conflict 'S114L4A5_13650';
%  CmdFileName={'S119L4A5_14391'};

%Fastigial & Vermis
%CmdFileName={'R167L1A3_20671';'R166L1A3_20371'}

for FileNb=1:length(CmdFileName);
    try
%% subject and procdir
if strcmp('R',CmdFileName{FileNb}(1))
   subject='Rigel';
   procdir = [directory,'processed',slash,'Rigel',slash];
elseif strcmp('S',CmdFileName{FileNb}(1))
   subject='Sixx';
   procdir = [directory,'processed',slash,'Sixx',slash];
elseif strcmp('H',CmdFileName{FileNb}(1))
   subject='Hilda';
   procdir = [directory,'processed',slash,'Hilda',slash];
end

procdirlisting=dir(procdir);
procdirfileNames={procdirlisting.name};
loadfile=procdirfileNames(~cellfun('isempty',regexpi(procdirfileNames,CmdFileName{FileNb},'match')));
% if two versions (REX and Spike2), chose Spike2
if length(loadfile)>1
    loadfile=loadfile(~cellfun('isempty',regexpi(loadfile,'Sp2','match')));
    if isempty(loadfile) %actually it's an old file that sneaked in
        loadfile=procdirfileNames(~cellfun('isempty',regexpi(procdirfileNames,CmdFileName{FileNb},'match')));
        loadfile=loadfile(~cellfun('isempty',regexpi(loadfile,'REX','match')));
    end
end

%% get countermanding session results 
[mssrt,inhibfun,ccssd,nccssd,ssdvalues,tachomc,tachowidth,sacdelay,rewtimes]=findssrt(loadfile{:}, 0);

%% align rasters
% set presets
tasktype='gapstop';
[~, trialdirs] = data_info(loadfile{:}, 1, 1); %reload file: yes (shouldn't happen, though), skip unprocessed files: yes

%alignments=1:3;

% sac vs stop
% firstalign=6; 
% secondalign=8; 
% aligntype='failed_fast';
% triplot = 0; 

% tgt vs stop
% firstalign=7; 
% secondalign=8; 
% aligntype='correct_slow';
% triplot = 1; 

% ssd
firstalign=7; % as if align to target
secondalign=507; 
aligntype='ssd';
triplot = 0; 

includebad=0;
spikechannel=1; %select appropriate cluster 
keepdir='compall'; %alldir %which sac directions
togrey=[];
singlerastplot=0;

if strcmp(aligntype,'ssd')
% if aligning to ssd, got to align NSS trials according to latency 
%allssds=unique([ccssd;nccssd]);

% canceled trials
ccssdval=unique(ccssd);
ctmatchlatidx=zeros(length(sacdelay),length(ccssdval));
for ssdval=1:length(ccssdval)
ctmatchlatidx(:,ssdval)=sacdelay>ccssdval(ssdval)+round(mssrt);
end
nullidx=sum(ctmatchlatidx,2)==0;
ctmatchlatidx(nullidx,1)=1;
% getting ssds for each NNS trial, taking the highest ssd. 
ctmatchlatidx=ccssdval(sum(ctmatchlatidx,2));
ctmatchlatidx(nullidx,1)=0;

% non-canceled trials
nccssdval=sort(unique(nccssd));
nctallmatchlatidx=zeros(length(sacdelay),length(nccssdval));
for ssdval=1:length(nccssdval)
nctallmatchlatidx(:,ssdval)=sacdelay>nccssdval(ssdval)+50 & sacdelay<nccssdval(ssdval)+round(mssrt);
end
% getting ssds for each NNS trial, taking the lowest ssd. 
nctmatchlatidx=zeros(size(nctallmatchlatidx,1),1);
for midx=1:size(nctallmatchlatidx,1)
    if ~isempty(find(nctallmatchlatidx(midx,:),1))
        nctmatchlatidx(midx)=nccssdval(find(nctallmatchlatidx(midx,:),1));
    end
end
end

%% use GUI-independent prealign
getaligndata = prealign(loadfile{:}(1:end-4), trialdirs, tasktype, firstalign,...
     secondalign,  includebad, spikechannel, keepdir,...
     togrey, singlerastplot, [ctmatchlatidx nctmatchlatidx]); % align data, don't plot rasters


%% plots
disp_cmd([loadfile{:}(1:end-4),'_Clus',num2str(spikechannel)],getaligndata,aligntype,triplot);

%     if strcmp(aligntype,'failed_fast') 
%             %     [p_cancellation,h_cancellation] = cmd_wilco_cancellation(rdd_filename,datalign);
%             disp_cmd([loadfile{:}(1:end-4),'_Clus',num2str(spikechannel)],getaligndata,aligntype,triplot);%0, 0: latmatch, no; triplot, no
%             %     disp_cmd(rdd_filename,datalign,1);
%     elseif strcmp(aligntype,'correct_slow') 
%             disp_cmd([loadfile{:}(1:end-4),'_Clus',num2str(spikechannel)],getaligndata,aligntype,triplot); % keep triplot off until fixed
%     elseif strcmp(aligntype,'ssd') % may need task-specific analysis
% 
%     else
%         
%      end
    catch
        continue
    end
end
%fsigma, mstart, mstop,

