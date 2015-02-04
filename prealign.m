 function datalign = prealign(filename, trialdirs, tasktype, firstalign,...
     secondalign,  includebad, spikechannel, keepdir,...
     togrey, singlerastplot, option)
 
% based on rdd_rasters_sdf, made GUI independent
% the option argument passes some additional data

global directory slash;

if nargin < 11
    option = [];
end
            
alignsacnum=0;
alignseccodes=[];
alignlabel=[];
secalignlabel=[];
collapsecode=0;
optiondat=[];

[fixcode, fixoffcode, tgtcode, tgtoffcode, saccode, ...
    stopcode, rewcode, ~, errcode1, errcode2, ~, basecode] = taskfindecode(tasktype);

%% first code:
if  firstalign==6 % mainsacalign button
    ecodealign=saccode;
        alignlabel='sac';
elseif firstalign==7 % tgtshownalign button
    ecodealign=tgtcode;
        alignlabel='tgt';
elseif firstalign==4 % rewardnalign button
    ecodealign=rewcode;
        alignlabel='rew';
elseif firstalign==8 % stopsignalign button
    if ~strcmp(tasktype,'gapstop')
        %faire qqchose. Ou pas.
    else
        ecodealign=stopcode;
            alignlabel='stop';
    end
elseif firstalign==9 % other sac align
    %made for corrective saccades in particular
    ecodealign=saccode;
    alignsacnum=1; %that is n-th saccade following the alignment code, which for now will be the main saccade.
        alignlabel='corsac';
else
    ecodealign=firstalign;
     if ecodealign==421
        alignlabel='touchbell';
    elseif ecodealign==742
        alignlabel='retarget';
    elseif ecodealign==507
        alignlabel='ssd';
    else
        alignlabel='ecode';
    end
end

%% second code: 
if secondalign==4
    secondcode=errcode1;
            secalignlabel='error1';
elseif secondalign==5
    secondcode=errcode2;
            secalignlabel='error2';
elseif secondalign==6
    secondcode=saccode;
            secalignlabel='sac';
elseif secondalign==7
    secondcode=tgtcode;
            secalignlabel='tgt';
elseif secondalign==8
    if ~strcmp(tasktype,'gapstop')
        %not good!
    else
        secondcode=stopcode;
            secalignlabel='stop';
        if firstalign==7 % tgtshownalign button
            if isnan(option)
                ecodealign=ecodealign(1); % no need to align stop trials to target, it will be done later
            elseif strcmp(option,'truealign')
                ecodealign=tgtcode(1); %except if it's really more convenient
                secondcode=tgtcode(2);
            end
        end
    end
elseif secondalign==9
    secondcode=saccode;
    alignsacnum=1;
            secalignlabel='corsac';
else
    secondcode=secondalign;
    if secondcode==507
            secalignlabel='ssd';
        if firstalign==7 % tgtshownalign button
            ecodealign(2)=ecodealign(1); %need to align to NSS trials twice
        end
    else   
            secalignlabel='none';
    end
end

%% Fusing task type and direction into ecode
if ecodealign<1000 % if only three numbers
    for i=1:length(trialdirs)
        aligncodes(i,:)=ecodealign*10+trialdirs(i);
    end
else
    aligncodes=ecodealign;
end
if logical(secondcode)
    if secondcode<1000
        for i=1:length(trialdirs)
            alignseccodes(i,:)=secondcode*10+trialdirs(i);
        end
    else
        alignseccodes=secondcode;
    end
end
if length(trialdirs)>1
    basecodes=[];
    for numbasecd=1:length(basecode)
        basecodes=[basecodes;(basecode(numbasecd)*ones(length(trialdirs),1)*10)+trialdirs];
    end
else
    basecodes=basecode;
end

if ~strcmp(keepdir,'alldir') && ~strcmp(keepdir,'compall');
  
    if strcmp(keepdir,'Horizontals') %remember error on initial gapstops
        codeidx=find(aligncodes-(floor(aligncodes./10).*10)==2 | aligncodes-(floor(aligncodes./10).*10)==6);
        aligncodes=aligncodes(codeidx);
        seccodeidx=find(alignseccodes-(floor(alignseccodes./10).*10)==2 | alignseccodes-(floor(alignseccodes./10).*10)==6);
        alignseccodes=alignseccodes(seccodeidx);
    elseif strcmp(keepdir,'Verticals')
        codeidx=find(aligncodes-(floor(aligncodes./10).*10)==0 | aligncodes-(floor(aligncodes./10).*10)==4);
        aligncodes=aligncodes(codeidx);
        seccodeidx=find(alignseccodes-(floor(alignseccodes./10).*10)==0 | alignseccodes-(floor(alignseccodes./10).*10)==4);
        alignseccodes=alignseccodes(seccodeidx);
    elseif strcmp(keepdir,'SU')
        codeidx=find(aligncodes-(floor(aligncodes./10).*10)==0);
        aligncodes=aligncodes(codeidx);
        seccodeidx=find(alignseccodes-(floor(alignseccodes./10).*10)==0);
        alignseccodes=alignseccodes(seccodeidx);
    elseif strcmp(keepdir,'UR')
        codeidx=find(aligncodes-(floor(aligncodes./10).*10)==1);
        aligncodes=aligncodes(codeidx);
        seccodeidx=find(alignseccodes-(floor(alignseccodes./10).*10)==1);
        alignseccodes=alignseccodes(seccodeidx);
    elseif strcmp(keepdir,'SR')
        codeidx=find(aligncodes-(floor(aligncodes./10).*10)==2);
        aligncodes=aligncodes(codeidx);
        seccodeidx=find(alignseccodes-(floor(alignseccodes./10).*10)==2);
        alignseccodes=alignseccodes(seccodeidx);
    elseif strcmp(keepdir,'BR')
        codeidx=find(aligncodes-(floor(aligncodes./10).*10)==3);
        aligncodes=aligncodes(codeidx);
        seccodeidx=find(alignseccodes-(floor(alignseccodes./10).*10)==3);
        alignseccodes=alignseccodes(seccodeidx);
    elseif strcmp(keepdir,'SD')
        codeidx=find(aligncodes-(floor(aligncodes./10).*10)==4);
        aligncodes=aligncodes(codeidx);
        seccodeidx=find(alignseccodes-(floor(alignseccodes./10).*10)==4);
        alignseccodes=alignseccodes(seccodeidx);
    elseif strcmp(keepdir,'BL')
        codeidx=find(aligncodes-(floor(aligncodes./10).*10)==5);
        aligncodes=aligncodes(codeidx);
        seccodeidx=find(alignseccodes-(floor(alignseccodes./10).*10)==5);
        alignseccodes=alignseccodes(seccodeidx);
    elseif strcmp(keepdir,'SL')
        codeidx=find(aligncodes-(floor(aligncodes./10).*10)==6);
        aligncodes=aligncodes(codeidx);
        seccodeidx=find(alignseccodes-(floor(alignseccodes./10).*10)==6);
        alignseccodes=alignseccodes(seccodeidx);
    elseif strcmp(keepdir,'UL')
        codeidx=find(aligncodes-(floor(aligncodes./10).*10)==7);
        aligncodes=aligncodes(codeidx);
        seccodeidx=find(alignseccodes-(floor(alignseccodes./10).*10)==7);
        alignseccodes=alignseccodes(seccodeidx);
    end
    
elseif strcmp(keepdir,'compall');
    collapsecode=1;
    aligncodes=aligncodes'; % previously: ecodealign,  so that when
    alignseccodes= alignseccodes'; %secondcode;
elseif strcmp(tasktype,'base2rem50')
    collapsecode=0; % we want to plots, one for each type of code
    aligncodes=aligncodes';
    alignseccodes= alignseccodes';
else
    disp('Selected option: all directions'); %correspond to 'selecalldir' tag
end

%% Grey area in raster
greycodes=[];
failcode=[errcode1 errcode2];

if strcmp(tasktype,'gapstop') %otherwise CAT arguments dimensions are not consistent below
    saccode=[saccode saccode];
    stopcode=[stopcode stopcode];
    rewcode=[rewcode rewcode];
    failcode=[errcode1 errcode1 errcode2 errcode2];
end

conditions =[tgtcode tgtoffcode;saccode saccode;fixcode fixoffcode;rewcode rewcode;failcode];

if logical(sum(togrey))
    greycodes=conditions(togrey,:); %selecting out the codes
end

%% Task-specific instructions
ol_instruct='directions'; %default mode
if strcmp(tasktype,'optiloc')
    return
%     ol_instructs=get(findobj('Tag','optiloc_popup'),'String');
%     ol_instruct=ol_instructs{get(findobj('Tag','optiloc_popup'),'Value')};
end

%% aligning data 

% default nonecodes. Potential conflict resolved in rdd_rasters
nonecodes=[17385 16386];

% variable to save aligned data
datalign=struct('dir',{},'rasters',{},'trials',{},'trigtosac',{},'sactotrig',{},...
    'trigtovis',{},'vistotrig',{},'alignidx',{},'eyeh',{},'eyev',{},'eyevel',{},...
    'amplitudes',{},'peakvels',{},'peakaccs',{},'allgreyareas',{},'stats',{},...
    'alignlabel',{},'savealignname',{},'bad',{},'extras',{});

if  singlerastplot || aligncodes(1)==1030 || aligncodes(1)== 17385
    datalign(1).alignlabel=alignlabel; %only one array
else
    for numlab=1:size(aligncodes,1)+size(alignseccodes,1)
        datalign(numlab).alignlabel =alignlabel;
    end
end

if sum(alignseccodes)    
    for numlab=size(aligncodes,1)+1:size(aligncodes,1)+size(alignseccodes,1)
        datalign(numlab).alignlabel =secalignlabel;
    end
end

%% formatting aligncodes

allaligncodes=[];

if ~sum(alignseccodes) %only one align code
    numcodes=size(aligncodes,1);
    allaligncodes=aligncodes;
    rotaterow=0;
else
    if collapsecode
        numcodes=size(aligncodes,1)+size(alignseccodes,1); %if collapsed ecodes
    else
        numcodes=2*max(size(aligncodes,1),size(alignseccodes,1)); %not collapsed together
    end
    if length(aligncodes)==length(alignseccodes)
        allaligncodes=[aligncodes;alignseccodes];
        rotaterow=0;
    else %unequal length of alignment codes. Making them equal here
        allaligncodes=1001*ones(numcodes,2); %first making a matrix 1001 to fill up the future "voids"
        if size(aligncodes,1)>size(alignseccodes,1)
            allaligncodes(1:size(aligncodes,1),1)=aligncodes;
            allaligncodes(size(aligncodes,1)+1:end,1)=alignseccodes*ones(size(aligncodes,1),1);
            allaligncodes(size(aligncodes,1)+1:end,2)=basecodes;
            rotaterow=fliplr(allaligncodes(size(aligncodes,1)+1:end,:));
        elseif size(aligncodes,1)<size(alignseccodes,1)
            allaligncodes(1:size(alignseccodes,1),1)=alignseccodes;
            allaligncodes(size(alignseccodes,1)+1:end,1)=aligncodes*ones(size(alignseccodes,1),1);
            allaligncodes(size(alignseccodes,1)+1:end,2)=basecodes;
            rotaterow=fliplr(allaligncodes(size(alignseccodes,1)+1:end,:));
        elseif size(aligncodes,2)>size(alignseccodes,2)
            allaligncodes(1:size(aligncodes,1),1:size(aligncodes,2))=aligncodes;
            allaligncodes(size(aligncodes,1)+1:end,1:size(alignseccodes,2))=alignseccodes*ones(size(aligncodes,1),1);
            allaligncodes(size(aligncodes,1)+1:end,size(alignseccodes,2)+1:end)=NaN;
            rotaterow=0;
        elseif size(aligncodes,2)<size(alignseccodes,2)
            allaligncodes(1:size(alignseccodes,1),1:size(alignseccodes,2))=alignseccodes;
            allaligncodes(size(alignseccodes,1)+1:end,1:size(aligncodes,2))=aligncodes*ones(size(alignseccodes,1),1);
            allaligncodes(size(alignseccodes,1)+1:end,size(aligncodes,2)+1:end)=NaN;
            rotaterow=0;
        end
    end
end

if strcmp(tasktype,'optiloc')
    if strcmp(ol_instruct,'directions') 
        %default, nothing to change
    elseif strcmp(ol_instruct,'amplitudes') && singlerastplot
        singlerastplot=0;
    elseif strcmp(ol_instruct,'directions and amplitudes') && singlerastplot
        singlerastplot=0;
    end
end

% align trials
for cnc=1:numcodes
    aligntype=datalign(cnc).alignlabel;
    adjconditions=conditions;
    if strcmp(tasktype,'gapstop')
        if strcmp(aligntype,'stop')
            includebad=1; %we want to include non-cancelled
            d_increment=size([aligncodes alignseccodes],1);%make room for additional "non-cancel" data
            numplots=numcodes+d_increment;
            if ~isempty(extras)
                optiondat=extras;
            elseif~isempty(option)
                optiondat=option;
            end
        elseif strcmp(aligntype,'ssd')
            includebad=1; %we want to compare cancelled with non-cancelled
            d_increment=size(alignseccodes,1);%make room for additional "non-cancel" data
            numplots=numcodes+d_increment;
        elseif strcmp(aligntype,'tgt') && strcmp(secalignlabel,'ssd')
                optiondat=option(:,cnc);
        elseif strcmp(aligntype,'rew')
            allaligncodes=flipud(allaligncodes); %the 1030 reward code should come first. It is not here for legacy reasons 
        end            
    elseif strcmp(tasktype,'base2rem50')
        adjconditions=[conditions(cnc,:);conditions(cnc+numcodes,:);conditions(cnc+2*numcodes,:)];
        numplots=numcodes;
    elseif logical(sum(strfind(ol_instruct,'amplitudes')))
        numplots=numcodes*3;
    else
       % includebad=0;
        numplots=numcodes;
    end 
    % Should simplify to structures rasters, indices, eyevalues, extras ...
    [rasters,aidx, trialidx, trigtosacs, sactotrigs, trigtovis, vistotrigs, eyeh,eyev,eyevel,...
        amplitudes,peakvels,peakaccs,allgreyareas,badidx,ssd,dirs,extras] = alignrasters( filename, tasktype, spikechannel, ...
        allaligncodes(cnc,:), nonecodes, includebad, alignsacnum, aligntype, collapsecode, adjconditions, firstalign,...
        optiondat);
    
    if isempty( rasters )
        disp( 'No raster could be generated (rex_rasters_trialtype returned empty raster)' );
        continue;
    elseif strcmp(aligntype,'stop') || strcmp(aligntype,'ssd')
        canceledtrials=~badidx';
        datalign(cnc).alignlabel='stop_cancel';
        datalign(cnc).rasters=rasters(canceledtrials,:);
        datalign(cnc).alignidx=aidx;
        datalign(cnc).dir=dirs(canceledtrials);
        datalign(cnc).trials=trialidx(canceledtrials);
        datalign(cnc).trigtosac=trigtosacs(canceledtrials);
        datalign(cnc).sactotrig=sactotrigs(canceledtrials);
        datalign(cnc).trigtovis=trigtovis(canceledtrials);
        datalign(cnc).vistotrig=vistotrigs(canceledtrials);
        datalign(cnc).eyeh=eyeh(canceledtrials,:);
        datalign(cnc).eyev=eyev(canceledtrials,:);
        datalign(cnc).eyevel=eyevel(canceledtrials,:);
        datalign(cnc).allgreyareas=allgreyareas(:,canceledtrials);
        datalign(cnc).amplitudes=amplitudes(canceledtrials);
        datalign(cnc).peakvels=peakvels(canceledtrials);
        datalign(cnc).peakaccs=peakaccs(canceledtrials);
        datalign(cnc).bad=badidx(canceledtrials);
        datalign(cnc).ssd=ssd(canceledtrials);
        
        canceledtrials=~canceledtrials;
        datalign(cnc+d_increment).alignlabel='stop_non_cancel';
        datalign(cnc+d_increment).rasters=rasters(canceledtrials,:);
        datalign(cnc+d_increment).alignidx=aidx;
        datalign(cnc+d_increment).dir=dirs(canceledtrials);
        datalign(cnc+d_increment).trials=trialidx(canceledtrials);
        datalign(cnc+d_increment).trigtosac=trigtosacs(canceledtrials);
        datalign(cnc+d_increment).sactotrig=sactotrigs(canceledtrials);
        datalign(cnc+d_increment).trigtovis=trigtovis(canceledtrials);
        datalign(cnc+d_increment).vistotrig=vistotrigs(canceledtrials);
        datalign(cnc+d_increment).eyeh=eyeh(canceledtrials,:);
        datalign(cnc+d_increment).eyev=eyev(canceledtrials,:);
        datalign(cnc+d_increment).eyevel=eyevel(canceledtrials,:);
        datalign(cnc+d_increment).allgreyareas=allgreyareas(:,canceledtrials);
        datalign(cnc+d_increment).amplitudes=amplitudes(canceledtrials);
        datalign(cnc+d_increment).peakvels=peakvels(canceledtrials);
        datalign(cnc+d_increment).peakaccs=peakaccs(canceledtrials);
        datalign(cnc+d_increment).bad=badidx(canceledtrials);
        datalign(cnc+d_increment).ssd=ssd(canceledtrials);
        %             datalign(cnc+d_increment).condtimes=condtimes(canceledtrials);
        elseif strcmp(tasktype,'optiloc') && logical(sum(strfind(ol_instruct,'amplitudes')))
        % compare amp distrib with expected distrib, typically [4,12,20]
        if ~sum(hist(abs(amplitudes),[4,12,20])==hist(abs(amplitudes),3))==3 %case when amps are not dixtributed as expected
            [~, apmbounds]=hist(abs(amplitudes),3);
            disp('unequal amp distrib in rdd_rasters_sdf line 515');
            pause;
        else
            shortampslim=4;
            medampslim=12;
            longampslim=20;
        end
        %allamps=(sort(abs(amplitudes)));
        
        shortamps=abs(amplitudes)<shortampslim; %(abs(amplitudes)<allamps(apmdistrib(1)))';
        datalign(cnc).alignlabel=[alignlabel,'4dg'];
        datalign(cnc).rasters=rasters(shortamps,:);
        datalign(cnc).alignidx=aidx;
        datalign(cnc).trials=trialidx(shortamps);
        datalign(cnc).trigtosac=trigtosacs(shortamps);
        datalign(cnc).sactotrig=sactotrigs(shortamps);
        datalign(cnc).trigtovis=trigtovis(shortamps);
        datalign(cnc).vistotrig=vistotrigs(shortamps);
        datalign(cnc).eyeh=eyeh(shortamps,:);
        datalign(cnc).eyev=eyev(shortamps,:);
        datalign(cnc).eyevel=eyevel(shortamps,:);
        datalign(cnc).allgreyareas=allgreyareas(:,shortamps);
        datalign(cnc).amplitudes=amplitudes(shortamps);
        datalign(cnc).peakvels=peakvels(shortamps);
        datalign(cnc).peakaccs=peakaccs(shortamps);
        datalign(cnc).bad=badidx(shortamps);
        
        medamps=abs(amplitudes)<medampslim;
        datalign(cnc+numcodes).alignlabel=[alignlabel,'12dg'];
        datalign(cnc+numcodes).rasters=rasters(medamps,:);
        datalign(cnc+numcodes).alignidx=aidx;
        datalign(cnc+numcodes).trials=trialidx(medamps);
        datalign(cnc+numcodes).trigtosac=trigtosacs(medamps);
        datalign(cnc+numcodes).sactotrig=sactotrigs(medamps);
        datalign(cnc+numcodes).trigtovis=trigtovis(medamps);
        datalign(cnc+numcodes).vistotrig=vistotrigs(medamps);
        datalign(cnc+numcodes).eyeh=eyeh(medamps,:);
        datalign(cnc+numcodes).eyev=eyev(medamps,:);
        datalign(cnc+numcodes).eyevel=eyevel(medamps,:);
        datalign(cnc+numcodes).allgreyareas=allgreyareas(:,medamps);
        datalign(cnc+numcodes).amplitudes=amplitudes(medamps);
        datalign(cnc+numcodes).peakvels=peakvels(medamps);
        datalign(cnc+numcodes).peakaccs=peakaccs(medamps);
        datalign(cnc+numcodes).bad=badidx(medamps);

        longamps=abs(amplitudes)<longampslim;
        datalign(cnc+2*numcodes).alignlabel=[alignlabel,'20dg'];
        datalign(cnc+2*numcodes).rasters=rasters(longamps,:);
        datalign(cnc+2*numcodes).alignidx=aidx;
        datalign(cnc+2*numcodes).trials=trialidx(longamps);
        datalign(cnc+2*numcodes).trigtosac=trigtosacs(longamps);
        datalign(cnc+2*numcodes).sactotrig=sactotrigs(longamps);
        datalign(cnc+2*numcodes).trigtovis=trigtovis(longamps);
        datalign(cnc+2*numcodes).vistotrig=vistotrigs(longamps);
        datalign(cnc+2*numcodes).eyeh=eyeh(longamps,:);
        datalign(cnc+2*numcodes).eyev=eyev(longamps,:);
        datalign(cnc+2*numcodes).eyevel=eyevel(longamps,:);
        datalign(cnc+2*numcodes).allgreyareas=allgreyareas(:,longamps);
        datalign(cnc+2*numcodes).amplitudes=amplitudes(longamps);
        datalign(cnc+2*numcodes).peakvels=peakvels(longamps);
        datalign(cnc+2*numcodes).peakaccs=peakaccs(longamps);
        datalign(cnc+2*numcodes).bad=badidx(longamps);
        
    else
        datalign(cnc).rasters=rasters;
        datalign(cnc).alignidx=aidx;
        datalign(cnc).dir=dirs;
        datalign(cnc).trials=trialidx;
        datalign(cnc).trigtosac=trigtosacs;
        datalign(cnc).sactotrig=sactotrigs;
        datalign(cnc).trigtovis=trigtovis;
        datalign(cnc).vistotrig=vistotrigs;
        datalign(cnc).eyeh=eyeh;
        datalign(cnc).eyev=eyev;
        datalign(cnc).eyevel=eyevel;
        datalign(cnc).allgreyareas=allgreyareas;
        datalign(cnc).amplitudes=amplitudes;
        datalign(cnc).peakvels=peakvels;
        datalign(cnc).peakaccs=peakaccs;
        datalign(cnc).bad=badidx;
        %             datalign(cnc).condtimes=condtimes;
    end
    
end


%% save name
if strcmp(tasktype,'optiloc')
    parsename=unique(cellfun(@(x) x(1:(regexp(x,'\d+')-1)), unique({datalign(~cellfun('isempty',{datalign.alignlabel})).alignlabel}), 'UniformOutput', false)); %e.g., {sac}, instead of {'sac12dg','sac20dg','sac4dg'}
else
    parsename=unique({datalign.alignlabel});
end

datalign(1).savealignname = cat( 2, directory, 'processed',slash, 'aligned',slash, filename, '_', cell2mat(parsename),'_c',num2str(spikechannel));

end
