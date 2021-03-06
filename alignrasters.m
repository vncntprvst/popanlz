function [alignedrasters, alignindex, trialindex, alltrigtosac, ...
    allsactotrig, alltrigtovis, allvistotrig,eyehoriz, eyevert, ....
    eyevelocity, sacspecs, allonoffcodetime,badidx,allssd,alldir,extras] = ...
    alignrasters( name, tasktype, selclus, aligntocode, noneofcodes,...
    allowbadtrials, alignsacnum, aligntype, collapse, conditions,...
    firstalign, option)

% based on rdd_rasters, but GUI independent

global rexnumtrials rexloadedname;

if strcmp(tasktype,'gapstop') || strcmp(tasktype,'base2rem50')
    multicodetask=1;
    if strcmp(aligntype,'tgt') && ~isempty(option)
        addshift = 1;
        shift=option; %that's the ssd, to realign from tgt to ss
        %     elseif strcmp(aligntype,'ssd')
        %          addshift = 1;
    else
        addshift = 0;
    end
else
    multicodetask=0;
end

if strcmp(aligntype,'stop') % get ssrt
    try
        mssrt=findssrt(name);
    catch
        disp('failed finding ssrt');
    end
end

%% variables init
bla=[];
[~, ~, tgtcode, tgtoffcode] = taskfindecode(tasktype);
alignedrasters=[];
sphisto=[];
alignindex=[];
trialindex=[];
trigtosac=[];
sactotrig=[];
trigtovis=[];
vistotrig=[];
eyehoriz=[];
eyevert=[];
eyevelocity=[];
trialnumbers=[];
%oldallonoffcodetime=[];

alignto = aligntocode;
sacamp = NaN;
sacpeakpeakvel = NaN;
sacpeakacc = NaN;
extras=[];

sbad = '';
if ~allowbadtrials
    sbad = 'Bad trials skipped: ';
end;

% Variables that will be incremented or appended as matching trials are
% collected.

alignmentfound = 0;
nummatch = 0;
alignindexlist = [];
rasters = [];
eyeh = [];
eyev = [];
eyevel = [];
eyehoriz = [];
eyevert = [];
eyevelocity = [];
allonoffcodetime = [];
onoffcodetime=[];
alltrigtosac=[];
allsactotrig=[];
alltrigtovis=[];
allvistotrig=[];
badidx=[];
allssd=[];
alldir=[];
% allcondtime = [];
sacspecs=struct('amplitudes',[],'peakvels',[],'peakaccs',[],'latency',[]);

cutcond=0;

%% Which Cluster?
%%%%%%%%%%%%%%%%%
% selclus = str2double(get(findobj('Tag','whichclus'),'String')); %
% unnecessary now: spikechannel and selclus are now one and the same

%  Loop through all of the trials using rex_first_trial and rex_next_trial.
%  See if each trial has the right codes, and try to align the spike data
%  to one of the alignment codes.

d = rex_first_trial( name, rexnumtrials, allowbadtrials );
islast = (d == 0);
while ~islast
    
    % For the current trial given by 'd', a call to rex_trial
    % gives us the codes, their times, the spike data, the
    % sampling rate, and the horizontal and vertical eye traces.  It also
    % gives start_time relative to the start of the whole file, and a
    % badtrial flag, which is irrelevant since we already know if this will
    % be a valid trial or not (because of the 2nd parameter in
    % rex_first_trial and rex_next_trial).
    
    %[ecodeout, etimeout, spkchan, spk, arate, h, v, start_time, badtrial ] = rex_trial(name, d );
    try %conditional breakpoint: logical(sum(ismember(ecodeout,alignto)))
        [ecodeout, etimeout, spkchan, spk, arate, h, v, start_time, isbadtrial, curtrialsacInfo] =...
            rdd_rex_trial(name, d, selclus);%, rdt_includeaborted);
        curdir=ecodeout(2)-floor(ecodeout(2)/10)*10;
%         if isempty(rexloadedname)
            rexloadedname=name;
%         end
    catch
        disp('call to rdd_rex_trial in alignrasters failed');
        [ecodeout, etimeout, spkchan, spk, arate, h, v, start_time, isbadtrial, curtrialsacInfo] = deal([]);
        curdir=NaN;
    end
    %if ~isbadtrial
    if logical(sum((ecodeout==2222)))  % had a weird case of a trial with ecode 2222. Don't know what that was. See file S110L4A5_12951 or S125L4A6_13990
        ecodeout(find(ecodeout==1035))=17385; % replace 1035 code by error code, to false positive on 1030
        isbadtrial=1;
    end
    %end
    
    
    %     if strcmp(tasktype,'gapstop') || strcmp(tasktype,'base2rem50')
    %         multicode=1;
    %     else
    %         multicode=0;
    %     end
    %
    
    if isempty(h) || isempty(ecodeout)
        cond_disp( 'Something wrong with trial, no data.' );
    else
        if  collapse || multicodetask%for collapsed alignements
            %             || multicode
            anyof = has_any_of( ecodeout, alignto );
            allof = 1;
        else
            allof = has_all_of( ecodeout, alignto );
            anyof=1;
        end
        if logical(sum(find(noneofcodes==alignto(1)))) ... %in case the purpose IS to align to a noneof code
                || allowbadtrials % or if we want the bad trials too
            noneof = 1;
            isbadtrial=~has_none_of(ecodeout, noneofcodes); %make sure trials with noneofcodes other than 17385 are tagged as bad
        else
            noneof = has_none_of( ecodeout, noneofcodes );
        end
        
        %  If these are all true, we have found a trial matching the
        %  requested codes.  Now check for alignment, which might be a
        %  whole list of possible candidates.  This actually makes a list
        %  (falign) of the ecode indices where there's a match, which is
        %  probably unneccessary, since only the first is used.
        
        if allof & noneof & anyof
            falign = [];
            for i = 1:length( alignto )
                fnext = find( ecodeout == alignto(i) );
                if ~isempty( fnext )
                    if isempty( falign )
                        falign = fnext(1);
                    else
                        falign = [falign;fnext(1)];
                    end
                end
            end
            
            if isempty( falign )
                s = sprintf( 'In rdd_rasters, trial %d has a matching base code (%d?), but does not contain any alignment code requested.', d, ecodeout(2) );
                cond_disp( s );
            else
                %% getting align times
                % We found one or more alignments, so get the actual time of the
                % first one. (we only want to align on the first one: make sure it's the right one passed as argument)
                
                alignmentfound = ecodeout( falign(1) );
                aligntime = etimeout( falign( 1 ) ) * (arate / 1000);
                
               if strcmp(aligntype,'error2') && ~(firstalign==4 || firstalign==5) %then purpose is not to align to error code
                   errortrial_align='yes'; 
               else
                   errortrial_align=[]; 
               end
                
                %% Adjust aligntime according to alignment type.
                %If it's a saccade align code, replace aligntime with the actual sac start (with new sac detection method:
                
                if  strcmp(aligntype,'sac') || strcmp(aligntype,'corsac') ...
                        || strcmp(aligntype,'error2')% mainsacalign button OR corrective saccade
                    ampsacofint=[];
%                     ampsacofint=zeros(1,size(curtrialsacInfo,2));
                    nwsacstart=cat(1,curtrialsacInfo.starttime);
                    if strcmp(tasktype,'tokens')
                        if strcmp(aligntype,'error2')
                            sacofint=nwsacstart>etimeout(falign(1)-2); % finding saccades occuring
                            % between the last token and
                            %the detected saccade (to the wrong target)
                        else
                            sacofint=nwsacstart>etimeout(falign(1))-40;  % the token task is special
                            % in that we do not detect the saccade itself, but the eye leaving
                            % the fixation window. The small delay (40ms) reflects that
                        end
                    else
                        sacofint=nwsacstart>etimeout(falign(1)-1); %considering all saccades occuring after the ecode
                        %preceding the saccade ecode, which is often erroneous
                    end
                    if strcmp(tasktype,'st_saccades')
                        if errortrial_align & firstalign==7
                            errortrial_align='tgt';
                        end
                        if nwsacstart(find(sacofint,1))-etimeout(7)<600 %there was a saccade during the delay period that went unnoticed by Rex saccade detection
%                             disp('check that code line in popanlz/alignraters:207')
%                             pause
                            if strcmp(aligntype,'error2') %that's fine, early saccade is supposed to happen 
                            elseif abs(curtrialsacInfo(find(sacofint,1)).amplitude)>2 % saving trials only if that delay-period 
                                                            %saccade was under 2 degrees, which happens pretty frequently with Sixx.
                                sacofint=0;
                                alignmentfound = 0;
                            end
                        end
                    end
%                     ampsacofint(find(sacofint,1):length(sacofint))=arrayfun(@(x) abs(getfield(curtrialsacInfo, {x}, 'amplitude')),find(sacofint,1):length(sacofint));
                    for k=find(sacofint,1):length(sacofint) %yes, this code is ugly.
                        ampsacofint(1,k)=abs(getfield(curtrialsacInfo, {k}, 'amplitude'));
                    end
                    %start time of first saccade greater than 2.5 degrees (typical
                    %restriction window) after relevant ecode (ecodesacstart-1)
                    if ~logical(sum(ampsacofint))
                        alignmentfound = 0;
                    elseif logical(sum(ampsacofint>2.5))
                        if strcmp(aligntype,'sac') || strcmp(aligntype,'error2')
                            if strcmp(errortrial_align,'tgt')
                                aligntime = etimeout(ismember(floor(ecodeout./10),tgtcode));
                            else
                                aligntime=getfield(curtrialsacInfo, {find(ampsacofint>2.5,1)}, 'starttime');
                            end
                            sacamp=getfield(curtrialsacInfo, {find(ampsacofint>2.5,1)}, 'amplitude');
                            sacpeakpeakvel=getfield(curtrialsacInfo, {find(ampsacofint>2.5,1)}, 'peakVelocity');
                            sacpeakacc=getfield(curtrialsacInfo, {find(ampsacofint>2.5,1)}, 'peakAcceleration');
                        elseif strcmp(aligntype,'corsac')  &&  find(ampsacofint>2.5,1)+1<=length(ampsacofint) % If we are looking for the n-th saccade after the main one
                            nextgoodsac=find(ampsacofint>2.5,1)+1;
                            aligntime=getfield(curtrialsacInfo, {nextgoodsac}, 'starttime');
                            if aligntime>etimeout(ecodeout==1030) || (getfield(curtrialsacInfo, {nextgoodsac}, 'starttime') - getfield(curtrialsacInfo, {nextgoodsac-1}, 'starttime')) > 250
                                alignmentfound = 0; % secondary saccade after reward or over 300 ms after main sac. Not considered corrective saccade
                            else
                                sacamp=getfield(curtrialsacInfo, {nextgoodsac}, 'amplitude');
                                sacpeakpeakvel=getfield(curtrialsacInfo, {nextgoodsac}, 'peakVelocity');
                                sacpeakacc=getfield(curtrialsacInfo, {nextgoodsac}, 'peakAcceleration');
                            end
                        else % no good main saccade or n-th saccade
                            alignmentfound = 0;
                        end
                    end
                end
                
                %% If it's a stop alignement code, two different cases: NCS or CS
                if strcmp(aligntype,'stop') && (sum(isnan(option)) || (firstalign==4 && logical(sum(option))))
                    if find(ecodeout==1503)
                        ncecode=10;
                    else
                        ncecode=9;
                    end
                    if isbadtrial % non-canceled stop trial
                        if firstalign==6
                            % for non-canceled stop trial, align to saccade
                            % initiation
                            if ecodeout(ncecode)==17385 || ecodeout(ncecode)==16386
                                nwsacstart=[curtrialsacInfo.starttime];
                                sacofint=nwsacstart>etimeout(7); %sacs that occur after event
                                sacofint=abs([curtrialsacInfo.amplitude])>2.5 & sacofint; % with a large enough amplitude
                                if sum(sacofint)
                                    aligntime= curtrialsacInfo(find(sacofint,1)).starttime;
                                else
                                    alignmentfound = 0;
                                end
                            end
                        elseif firstalign==7 % align to target
                            aligntime = etimeout(find(floor(ecodeout./10) == 487,1)) * (arate / 1000);
                        elseif firstalign==9 % align to corrective saccade after main (incorrect) saccade
                            nwsacstart=[curtrialsacInfo.starttime];
                            sacofint=nwsacstart>etimeout(7); %sacs that occur after event
                            sacofint=abs([curtrialsacInfo.amplitude])>2.5 & sacofint; % with a large enough amplitude
                            sacofint=find(sacofint, 2, 'first');
                            if length(sacofint)>1 && (curtrialsacInfo(sacofint(2)).starttime - curtrialsacInfo(sacofint(1)).starttime)<250
                                aligntime=curtrialsacInfo(sacofint(2)).starttime;
                                curtrialsacInfo(sacofint(2)).status='corsac';
                            else
                                alignmentfound = 0;
                            end
                         elseif firstalign==4 
                             %align to time of error code
%                              aligntime=etimeout(find(ecodeout==17385 | ecodeout==16386,1)) * (arate / 1000); 
                             % or
                             % align to expected reward time, based on Lm-NSS trials rew times
                             if ecodeout(ncecode)==17385 || ecodeout(ncecode)==16386
                                nwsacstart=[curtrialsacInfo.starttime];
                                sacofint=nwsacstart>etimeout(7); %sacs that occur after event
                                sacofint=abs([curtrialsacInfo.amplitude])>2.5 & sacofint; % with a large enough amplitude
                                if sum(sacofint)
                                    aligntime= curtrialsacInfo(find(sacofint,1)).starttime+round(mean(option)); 
                                else
                                    alignmentfound = 0;
                                end
                             end
                        end
                    else
                        % for successfully canceled stop trials, align to stop
                        % signal delay + stop signal reaction time, 
                        if firstalign==4 % unless the first alignement is to reward (then align to reward)
                            aligntime=etimeout(find(ecodeout == 1030,1)) * (arate / 1000);
                        elseif firstalign==7 % or target! then align to target ... 
                            aligntime=etimeout(find(floor(ecodeout./10) == 487,1)) * (arate / 1000);
                        else
                            if ~isnan(mssrt)
                                aligntime=etimeout( falign( 1 ) ) * (arate / 1000)+round(mssrt);
                            else
                                alignmentfound=0;
                            end
                        end
                    end
                end
                
                %% now get the time of "grey area" ecodes
                % used to get only times for selected conditions. Now get
                % all, and just displaying the requested ones in
                % rdd_rasters_sdf. VP 7/14/2012
                
                greytypes={'cue';'eyemvt';'fix';'rew';'fail'};
                failedsac=0;
                %caveat: some conditions may be 4 or 5 digits long,
                %such as user defined codes such as TOKSWCD (1501)
                if ~cutcond
                    conditions(conditions>=1000)=floor(conditions(conditions>=1000)/10); %cut last digit off of them
                    cutcond=1;
                end
                
                %                 greytypes=(greytypes(selectedgrey));
                codepairnb=floor(size(conditions,2)/2);%there may be multiple code. See l. 421 as well as here (273)
                %              if logical(sum(greycodes))
                shortecodout=floor(ecodeout./10);
                %                     if size(greycodes,1)>1 %more than one row: means several checkboxes are selected
                for i=1:size(conditions,1)
                    goodsacnum=0;
                    if strcmp(greytypes(i),'eyemvt') %adjust times to real saccade times
                        % find which saccade is the "good" one (if any) in this trial
                        try 
                            if length(find(~cellfun(@isempty,{curtrialsacInfo.latency})))==1
                                goodsacnum=find(~cellfun(@isempty,{curtrialsacInfo.latency}));
                            elseif length(find(~cellfun(@isempty,{curtrialsacInfo.latency})))>1 && ... % toomany good sac
                                     ... %diff(find(~cellfun(@isempty,{curtrialsacInfo.latency})))==1 && ... % are consecutive
                                    abs(diff([curtrialsacInfo(~cellfun(@isempty,{curtrialsacInfo.latency})).latency]))<50 % under 50ms diffference
                                goodsacnum=find(~cellfun(@isempty,{curtrialsacInfo.latency}),1,'last');
                            elseif isempty(find(~cellfun(@isempty,{curtrialsacInfo.latency}), 1))
                                goodsacnum=0;
                            else
                                %problem. Got to fix REX_process
                             goodsacnum=find(~cellfun(@isempty,{curtrialsacInfo.latency}))   
                            end
                        catch
                            goodsacnum=0;
                        end
                        if ~logical(sum(goodsacnum)) && (~strcmp(aligntype,'stop') && ~strcmp(aligntype,'ssd') && ~strcmp(aligntype,'touchbell'))
%                             s = sprintf('alignrasters: cannot display grey area for trial %d because saccade cannot be found. Removing erroneous trial',d);
%                             disp(s);%ploteyesac(d,name)
                            %                             pause
                            failedsac=1;
                            alignmentfound = 0;
                        elseif strcmp(aligntype,'rew') && alignmentfound &&... % case of reward aligned where the "good" 
                                curtrialsacInfo(1,goodsacnum(1)).starttime > etimeout(ecodeout==alignmentfound) % saccade occurs after reward
                            failedsac=1;
                            alignmentfound = 0; %same punishment: out
                        end
                    end
                    for j=1:size(conditions,2)
                        try
                            fonoffcode(i,j) = find(shortecodout == conditions(i,j),1);
                        catch
                            fonoffcode(i,j) = NaN;
                        end %in multiple code tasks (eg, gapstop), may fail to find the code in the ecode list
                        try
                            onoffcodetime(i,j) = etimeout(fonoffcode(i,j)) * (arate / 1000);
                        catch
                            onoffcodetime(i,j) = NaN;
                        end
                    end
                    if logical(goodsacnum)
                        onoffcodetime(i,1)=getfield(curtrialsacInfo, {goodsacnum}, 'starttime');
                        onoffcodetime(i,1+codepairnb)=getfield(curtrialsacInfo, {goodsacnum}, 'endtime');
                    end
                    if i==4 && onoffcodetime(4,1)-onoffcodetime(1,1)>900 %Hilda's recordings
%                         onoffcodetime 
                    end
                end
                
                %% former code
                % else % if only one condition selected
                %                         goodsacnum=0;
                %                         if strcmp(greytypes,'eyemvt') %adjust times to real saccade times
                %                                 % find which saccade is the "good" one (if any) in this trial
                %                         try goodsacnum=find(~cellfun(@isempty,{curtrialsacInfo.latency}));   catch goodsacnum=0; end
                %                         end
                %                             for j=1:size(greycodes,2)
                %                                 try fonoffcode(j) = find( shortecodout == greycodes(j),1); catch fonoffcode(j) = NaN; end
                %                                 try onoffcodetime(j) = etimeout(fonoffcode(j)) * (arate / 1000); catch onoffcodetime(j) = NaN; end
                %                             end
                %                                if logical(goodsacnum)
                %                                     onoffcodetime(1)=getfield(curtrialsacInfo, {goodsacnum}, 'starttime');
                %                                     onoffcodetime(1+codepairnb)=getfield(curtrialsacInfo, {goodsacnum}, 'endtime');
                %                                 end
                %                     end
                %              end
                
                %% find condition times
                %                 codepairnb=floor(size(conditions,2)/2);%there may be multiple code. See l. 421 as well as here (273)
                %                 shortecodout=floor(ecodeout./10);
                %                 for i=1:size(conditions,1)
                %
                %                     for j=1:size(conditions,2)
                %                         try fcondcode(i,j) = find(shortecodout == conditions(i,j),1);  catch fcondcode(i,j) = NaN; end %in multiple code tasks (ie, gapstop), may fail to find the code in the ecode list
                %                         try condcodetime(i,j) = etimeout(fcondcode(i,j)) * (arate / 1000); catch condcodetime(i,j) = NaN; end
                %                     end
                %
                %                 end
                %
                %                 condmattime = {condcodetime};
                %
                %% filling up alignindexlist with align times
                %                if alignsacnum == 0
                
                % if there's a shift, may add it now
                if strcmp(tasktype,'gapstop') && strcmp(aligntype,'tgt') && addshift && ~failedsac
                    if ~shift(nummatch + 1)
                        bla=[bla; d];
                        alignmentfound = 0;
                        shfidx=true(size(shift));
                        shfidx(nummatch + 1)=0;
                        shift=shift(shfidx);
                    else
                        aligntime=aligntime+shift(nummatch + 1);
                    end
                end
                
                if alignmentfound
                    nummatch = nummatch + 1;
                    badidx(nummatch)=isbadtrial;
                    if strcmp(aligntype,'stop') || strcmp(aligntype,'ssd')
                        if ecodeout(9)==17385 || ecodeout(9)==16386
                            badidx(nummatch)=2; %non-cancelled trials
                        end
                        if ecodeout(8)==1503 %with photodiode state timestamp
                            allssd(nummatch,1)=etimeout(:,9)-etimeout(:,7)-3;
                            allssd(nummatch,2)=etimeout(:,8)-etimeout(:,7)-2;
                        else
                            allssd(nummatch,1)=etimeout(:,8)-etimeout(:,7)-2;
                        end
                        %                         if addshift %no need to add ssd shift, if it's not
                        %                                       %aligned to target
                        %                             aligntime=aligntime+allssd(nummatch);
                        %                         end
                    end
                    alignindexlist( nummatch ) = aligntime;
                    alldir( nummatch ) = curdir;
                    trialindex(nummatch)=d;
                    
                    % trigger times
                    if (tgtcode<=1000)
                        ecodeson=shortecodout;
                    else
                        ecodeson=ecodeout;
                    end
                    if (tgtoffcode<=1000)
                        ecodesoff=shortecodout;
                    else
                        ecodesoff=ecodeout;
                    end
                    visevents=[etimeout(ismember(ecodeson,tgtcode)), etimeout(ismember(ecodesoff,tgtoffcode))]-etimeout(1)-1;
                    
                    %collect "extras" 
                    if firstalign==4 && strcmp(aligntype,'rew') % keep time from saccade to reward
                        %get saccade time
                        nwsacstart=[curtrialsacInfo.starttime];
                        % Adjust timeout ecode if needed
                        sacecodenb=8; %Here assuming gapstop
                        sacofint=nwsacstart>etimeout(sacecodenb-1);
                        sacofint=abs([curtrialsacInfo.amplitude])>2.5 & sacofint; % with a large enough amplitude        
                        %start time of first saccade greater than 2.5 degrees (typical
                        %restriction window) after relevant ecode
                        if sum(sacofint)
%                             if aligntime-curtrialsacInfo(find(sacofint,1)).starttime < 0
%                                 curtrialsacInfo(find(sacofint,1)).starttime
%                             end
                            extras(nummatch,1)=aligntime-curtrialsacInfo(find(sacofint,1)).starttime;
                        else
                            disp('hum, why??') 
                        end
                    end
                    
                    
                    % recordings with trigger channel (sync' trigger, e.g. for spike2 recordings)
                    if find(ecodeout==1502) % Trigger code
                        %                         triggercode=1;
                        trigtosac=aligntime-etimeout(1)-1; %trigger code is 1ms before 1001
                        trigtovis=max(visevents(visevents<=trigtosac)); %the latest visual event occuring before alignment time
                        if find(ecodeout==1030)
                            sactotrig=etimeout(find(ecodeout==1502,1))+1-aligntime; %the second trigger channel is actually the start of the next trial
                            vistotrig=etimeout(find(ecodeout==1502,1))+1-max(visevents(visevents<=trigtosac)+etimeout(1)+1);
                        else
                            sactotrig=NaN;
                            vistotrig=NaN;
                        end
                        
                    else %older recordings without trigger code
                        %                         triggercode=0;
                        trigtosac=aligntime-etimeout(1)-1; %in case there is a trigger channel available in the SH recording
                        trigtovis=max(visevents(visevents<=trigtosac)); %the latest visual event occuring before alignment time
                        if find(ecodeout==1030) %good trial
                            sactotrig=etimeout(find(ecodeout==1030,1))+1-aligntime;%1ms between reward code and valve opening
                            vistotrig=etimeout(find(ecodeout==1030,1))+1-max(visevents(visevents<=trigtosac)+etimeout(1)+1);
                        else %wrong trial
                            sactotrig=NaN;
                            vistotrig=NaN;
                        end
                        if sactotrig<0
                            ecodeout;
                        end
                    end
                    
                end
                %%                 elseif alignsacnum < 0
                %                     cond_disp( 'In rdd_rasters, aligning to saccades BEFORE alignment codes has not been implemented yet. ');
                %                     alignmentfound = 0;
                %                 elseif alignsacnum > 0
                %                     [sstarts, sends] = rex_trial_saccade_times( name, d );
                %                     sacnumstart = find( sstarts > aligntime );
                %                     if length( sacnumstart ) < alignsacnum
                %                         alignmentfound = 0;
                %                     else
                %                         aligntime = sstarts( sacnumstart(alignsacnum) );
                %                         nummatch = nummatch + 1;
                %                         alignindexlist( nummatch ) = aligntime;
                %                     end;
                %                 end;
                
                %% Collecting spikes and stuff
                % If we found a place to align, either a code, or a code
                % followed by the alignsacnum-th saccade, then collect the
                % spikes for this trial in 'train', and then add to the
                % raster list of all spike trains so far, i.e. 'rasters'.
                % Though added to the rasters list, these trains are not
                % yet aligned.  That happens later.
                
                if alignmentfound
                    trialnumbers(nummatch)=d;
                    sacspecs(nummatch).amplitudes=sacamp;
                    sacspecs(nummatch).peakvels=sacpeakpeakvel;
                    sacspecs(nummatch).peakaccs=sacpeakacc;
                    % store latency as well 
                    if isempty([curtrialsacInfo.latency])
                        sacspecs(nummatch).latency=NaN; % likely CS trial
                    elseif ~exist('ampsacofint') 
                        sacspecs(nummatch).latency=[getfield(curtrialsacInfo, {find(arrayfun(@(x) ~isempty(x.latency), curtrialsacInfo) ,1)}, 'latency')];
                    else
                        try
                            sacspecs(nummatch).latency=[getfield(curtrialsacInfo, {find(ampsacofint>2.4,1)}, 'latency')];
                        catch
                            sacspecs(nummatch).latency=[];
                        end
                    end
                    train = [0];
                    if ~isempty( spk )
                        [train, last] = rex_spk2raster( spk, 1, max([aligntime+300, length(h)]) );
                        if isempty(train)
                            train=nan(1,length(h));
                            last=length(h);
                        end
                    end
                    if aligntime > length(train) %that's weird (added max([aligntime+300, length(h)] above)
%                         aligntime;
                    end
                    rasters = cat_variable_size_row(rasters, train);
                    
                    
                    %collect conditions (aka greycodes) times
                    if  ~(size(onoffcodetime,2)==2 && codepairnb==1)
                        onoffkeepcode=[find(~isnan(onoffcodetime(1,:)),1) find(~isnan(onoffcodetime(1,:)),1)+codepairnb];
                        if onoffkeepcode(end)>size(onoffcodetime,2)
                            onoffkeepcode=[size(onoffcodetime,2)-(codepairnb+1) size(onoffcodetime,2)-1];
                        end
                        onoffcodetime=onoffcodetime(:,onoffkeepcode);
                    end
                    allonoffcodetime=[allonoffcodetime {onoffcodetime}];
                    
                    %and collect trigger alignments
                    alltrigtosac=[alltrigtosac trigtosac];
                    allsactotrig=[allsactotrig sactotrig];
                    alltrigtovis=[alltrigtovis trigtovis];
                    allvistotrig=[allvistotrig vistotrig];
                    
                    
                    if length(h)<length(train)
                        s = sprintf( 'In rdd_rasters, the eye trace was shorter than the spike raster (%d < %d) for trial %d.  Padding with zeros.',...
                            length(h), length(train), d );
                        %cond_disp(s);
                        h = [h zeros(1, length(train)-length(h))];
                        v = [v zeros(1, length(train)-length(v))];
                    end;
                    
                    % Also collect eye movement traces for this trial, and
                    % add to the lists (eyeh and eyev).  Also do velocity.
                    
                    eyeh = cat_variable_size_row( eyeh, h );
                    eyev = cat_variable_size_row( eyev, v );
                    %                     dh = diff( h );
                    %                     dv = diff( v );
                    %                     velocity = sqrt( ( dh .* dh ) + ( dv .* dv ) );
                    [filth, filtv, filtvel]=cal_velacc(h,v);
                    eyevel = cat_variable_size_row( eyevel, filtvel);
                    
                    
                    %                     if logical(sum(greycodes))
                    %                         trialonofftime=zeros(1,length(h));
                    %                         if size(greycodes,1)>1 %more than one row
                    %                             for i=1:size(conditions,1)
                    %                                 onoffkeepcode=[onoffcodetime(i,find(~isnan(onoffcodetime(i,:)),1))...
                    %                                     onoffcodetime(i,find(~isnan(onoffcodetime(i,:)),1)+codepairnb)]; %keep the first pair of good codes
                    %                                trialonofftime(onoffkeepcode(1):onoffkeepcode(2))=i; % a bit ridiculous, but simpler than keeping indexes in this code
                    %                             end
                    %                         else
                    %                             onoffkeepcode=[onoffcodetime(1,find(~isnan(onoffcodetime(1,:)),1))...
                    %                             onoffcodetime(1,find(~isnan(onoffcodetime(1,:)),1)+codepairnb)]; %keep the first pair of good codes
                    %                             trialonofftime(onoffkeepcode(1):onoffkeepcode(2))=1;
                    %                         end
                    %
                    %                         allonoffcodetime=cat_variable_size_row(allonoffcodetime, trialonofftime);
                    %                     end
                    
                    % collect condition times
                    %                     allcondtime = cat_variable_size_row(allcondtime, condmattime);
                    
                    
                end
            end
        end
    end
    
    [d, islast] = rex_next_trial( name, d, allowbadtrials );
    
end


if isempty( rasters )
    cond_disp( 'rdd_rasters: Cannot generate rasters with the given codes, since no matching trials were found.' );
    alignedrasters = [];
    alignindex = 0;
    sphisto = [];
    return;
end

% We have rows of spike trains (rasters), and indices on which to align
% them (alignindexlist).  Now shift, or align, each of the rows so that the
% alignment time occurs at the same index in all rows.  See
% align_rows_on_indices() to see how this works.
alignindex = max( alignindexlist );

% alignindex is now the index (the column) in alignedrasters that is the
% column to which all rows are aligned.
[alignedrasters,shift] = align_rows_on_indices( rasters, alignindexlist );

%Do the same for the eye stuff.
eyehoriz = align_rows_on_indices( eyeh, alignindexlist );
eyevert = align_rows_on_indices( eyev, alignindexlist );
eyevelocity = align_rows_on_indices( eyevel, alignindexlist );

%add shift to grey areas times
for shifttm=1:size(shift,1)
    allonoffcodetime(shifttm)={cell2mat(allonoffcodetime(shifttm))+shift(shifttm)};
end
% cell2mat(allonoffcodetime(1))
% cell2mat(allonoffcodetime(1))+shift(1)
%
% foo = align_rows_on_indices(oldallonoffcodetime, alignindexlist );
%     find(oldallonoffcodetime(1,:),1)
%     find(foo(1,:),1)

end

% figure( 21 )
% subplot( 2, 1, 1 );
% imagesc( rasters );
% subplot( 2, 1, 2 );
% imagesc( alignedrasters );
% colormap( 1 - GRAY );
% alignindex
% alignindexlist

%same treatment for the on/off ecodes time
% if logical(sum(conditions))
%     allonofftime = align_rows_on_indices(allonoffcodetime, alignindexlist );
% else
%     allonofftime = [];
% end



% Done.  Everything is collected and aligned for all matching trials.

