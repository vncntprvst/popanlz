
filename='H125L6A2_17330'; % H125L6A2_17330 H125L6A2_16160
alignment='error2sac'; % 'sac';
cluster='1';

if strcmp(alignment,'error2sac') || strcmp(alignment,'sac')
    preblocks=6;
    postblocks=6;
end

% get overall crosscorrelogram, as well as correlation and coherence for
% x pre-event and y post-event blocks 
[cohrfreq,cohrmag,SFcorr,ntrials]=SFC(filename,alignment,cluster,preblocks,postblocks,1);

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
figure
subplot(2,1,1)
plot(cohrfreq{1,1},zcohrmag{1,1});
hold on
plot(cohrfreq{1,2},zcohrmag{1,2},'r');
ylim([0 0.2]);  xlim([3 50])
subplot(2,1,2)
plot(SFcorr{1});
hold on
plot(SFcorr{2},'r');
set(gca,'ylim',[-1 1],'xlim',[0 201],'xtick',1:100:201,'xticklabel',[-100 0 100])
