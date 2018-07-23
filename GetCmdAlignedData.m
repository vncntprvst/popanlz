function alignedData=GetCmdAlignedData(fileName, aligntype, globVars)
global rexnumtrials recDir
rexnumtrials = globVars{1};
recDir = globVars{2};
% recDir='E:\Data\Recordings\processed\Sixx\';
%% get countermanding session results
[mssrt,~,ccssd,nccssd,~,tachomc,~,sacdelay,tachowidth,...
    ~,~,~,trialdirs]=findssrt(fileName, 0);

tasktype='gapstop';

%alignments=1:3;
switch aligntype
    case 'saccade' %sac vs stop
        firstalign=6;
        secondalign=8;
        aligntype='failed_fast';
        option=NaN;
    case 'target'% tgt vs stop
        firstalign=7;
        secondalign=8;
        aligntype='correct_slow';
        %                         if plottype == 3;
        %                     option='truealign';
        %                 else
        option=NaN;
        %                 end
    case 'stopsignal' % stop signal
        firstalign=7; % as if align to target
        secondalign=507;
        aligntype='ssd';
        option=NaN;
        
        %     case 'corrsac' % align to corrective saccades in failed cancellation trial
        %                 firstalign=9; % corrective saccades for SST
        %                 option=NaN;
        
    case 'reward' % reward
        firstalign=4;
        secondalign=8;
        aligntype='rwd';
        option=NaN;
end

includebad=0;
spikechannel=1; %select appropriate cluster
keepdir='compall'; %alldir %which sac directions
togrey=[];
singlerastplot=0;

if strcmp(aligntype,'ssd')
    % if aligning to stop signal, got to align NSS trials according to latency
    %allssds=unique([ccssd;nccssd]);
    % canceled trials
    ccssdval=unique(ccssd);
    ctmatchlatidx=zeros(length(sacdelay.nsst),length(ccssdval));
    for ssdval=1:length(ccssdval)
        % Keeping NSS trials with sac latencies long enough
        % that they would have occured after a stop-signal
        ctmatchlatidx(:,ssdval)=sacdelay.nsst>ccssdval(ssdval)+round(mssrt);
    end
    nullidx=sum(ctmatchlatidx,2)==0;
    ctmatchlatidx(nullidx,1)=1;
    % getting ssds for each NNS trial, taking the highest ssd.
    ctmatchlatidx=ccssdval(sum(ctmatchlatidx,2));
    ctmatchlatidx(nullidx,1)=0;
    
    % non-canceled trials
    nccssdval=sort(unique(nccssd));
    nctallmatchlatidx=zeros(length(sacdelay.nsst),length(nccssdval));
    for ssdval=1:length(nccssdval) % keeping NSS trials in which
        % a saccade would have been
        % initiated even if a stop
        % signal had occurred, but
        % with saccade latencies
        % greater than the
        % stop-signal delay plus a
        % visual-response latency.
        % We take tachomc-tachowidth/2
        % rather than the arbitrary
        % 50ms from Hanes et al 98
        nctallmatchlatidx(:,ssdval)=sacdelay.nsst'>nccssdval(ssdval)+...
            (mean(tachomc)-tachowidth/2) & sacdelay.nsst'<nccssdval(ssdval)+round(mssrt);
    end
    % getting ssds for each NNS trial, taking the lowest ssd.
    nctmatchlatidx=zeros(size(nctallmatchlatidx,1),1);
    for midx=1:size(nctallmatchlatidx,1)
        if ~isempty(find(nctallmatchlatidx(midx,:),1))
            nctmatchlatidx(midx)=nccssdval(find(nctallmatchlatidx(midx,:),1));
        end
    end
    option=[ctmatchlatidx nctmatchlatidx];
    
    
    %% remove code commented below if all goes well
    %     % canceled trials
    %     ccssdval=unique(ccssd);
    %     ctmatchlatidx=zeros(length(sacdelay.nsst),length(ccssdval));
    %     for ssdval=1:length(ccssdval)
    %         ctmatchlatidx(:,ssdval)=sacdelay.nsst>ccssdval(ssdval)+round(mssrt);
    %     end
    %     nullidx=sum(ctmatchlatidx,2)==0;
    %     ctmatchlatidx(nullidx,1)=1;
    %     % getting ssds for each NNS trial, taking the highest ssd.
    %     ctmatchlatidx=ccssdval(sum(ctmatchlatidx,2));
    %     ctmatchlatidx(nullidx,1)=0;
    %
    %     % non-canceled trials
    %     nccssdval=sort(unique(nccssd));
    %     nctallmatchlatidx=zeros(length(sacdelay.nsst),length(nccssdval));
    %     for ssdval=1:length(nccssdval)
    %         nctallmatchlatidx(:,ssdval)=sacdelay.nsst>nccssdval(ssdval)+50 & sacdelay.nsst<nccssdval(ssdval)+round(mssrt);
    %     end
    %     % getting ssds for each NNS trial, taking the lowest ssd.
    %     nctmatchlatidx=zeros(size(nctallmatchlatidx,1),1);
    %     for midx=1:size(nctallmatchlatidx,1)
    %         if ~isempty(find(nctallmatchlatidx(midx,:),1))
    %             nctmatchlatidx(midx)=nccssdval(find(nctallmatchlatidx(midx,:),1));
    %         end
    %     end
end

%% use GUI-independent prealign
alignedData = prealign(fileName, trialdirs, tasktype, firstalign,...
    secondalign,  includebad, spikechannel, keepdir,...
    togrey, singlerastplot, option);