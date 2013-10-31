listing=dir('E:\Data\Recordings\processed\aligned');
session='H125L6A2';
alignment='error2sac'; % 'sac';
proc='all'; %REX or Sp2
cluster='all'; % cluster number, e.g. '_c1'
sesfilesidx=~cellfun('isempty',cellfun(@(x) strfind(x,session), {listing.name},'UniformOutput',false ));
noSHfiles=cellfun('isempty',cellfun(@(x) strfind(x,'2SH'), {listing.name},'UniformOutput',false ));
aligfilesidx=~cellfun('isempty',cellfun(@(x) strfind(x,alignment), {listing.name},'UniformOutput',false ));
if strcmp(proc,'all')
    procfilesidx=1;
else
    procfilesidx=~cellfun('isempty',cellfun(@(x) strfind(x,proc), {listing.name},'UniformOutput',false ));
end
if strcmp(proc,'all')
    clusfilesidx=1;
else
    clusfilesidx=~cellfun('isempty',cellfun(@(x) strfind(x,cluster), {listing.name},'UniformOutput',false )); 
end

FileList={listing(sesfilesidx & noSHfiles & aligfilesidx & procfilesidx & clusfilesidx).name};

% cross-correlation window
corrwind=100;

for filenm=1:length(FileList)
filename=FileList{filenm}; 
if strcmp(alignment,'error2sac') || strcmp(alignment,'sac')
    preblocks=6;
    postblocks=6;
end

% get overall crosscorrelogram, as well as correlation and coherence for
% x pre-event and y post-event blocks 
try
[cohrfreq,cohrmag,SFcorr,ntrials]=SFC(filename,alignment,cluster,corrwind,preblocks,postblocks,1);
catch
    continue;
end

numcomp=2;
zcohrmag=cell(2,numcomp);
zcohrmag_sem=cell(2,numcomp);
% Fisher z transform
for compnum=1:numcomp
zcohrmag{1,compnum}=.5.*log((1+cohrmag{1,compnum})./(1-cohrmag{1,compnum}));
zcohrmag_sem{1,compnum}=(1/sqrt(ntrials(compnum)-3))*1.96; %95 confidence interval
zcohrmag{2,compnum}=.5.*log((1+cohrmag{2,compnum})./(1-cohrmag{2,compnum}));
zcohrmag_sem{2,compnum}=(1/sqrt(ntrials(compnum)-3))*1.96; %95 confidence interval
end

%plots
figure('name',filename,'numbertitle','off')
subplot(2,1,1)
plot(cohrfreq{1,1},zcohrmag{1,1});
hold on
plot(cohrfreq{1,2},zcohrmag{1,2},'r');
ylim([0 0.2]);  xlim([3 50])    
title({'Coherence estimate'},'FontSize',20,'FontName','calibri');
xlabel({'Frequency (Hz)'},'FontSize',16,'FontName','calibri');
ylabel({'Magnitude'},'FontSize',16,'FontName','calibri');
subplot(2,1,2)
plot(SFcorr{1});
hold on
plot(SFcorr{2},'r');
set(gca,'ylim',[-1 1],'xlim',[0 2*corrwind+1],'xtick',1:corrwind:2*corrwind+1,'xticklabel',[-2*corrwind 0 2*corrwind])
title({'Spike-field correlation'},'FontSize',20,'FontName','calibri');
xlabel({'Time (ms)'},'FontSize',16,'FontName','calibri');
ylabel({'Magnitude'},'FontSize',16,'FontName','calibri');

end
