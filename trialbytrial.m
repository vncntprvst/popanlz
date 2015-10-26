function [behav, neur]=trialbytrial(gsdata,conn);

% global directory slash;
% if isempty(directory)
%     [directory,slash]=SetUserDir;
% end

[behav, neur]=deal(cell(1,1));

%% old code
% [peakcct, peaksdf,tbtdircor,tbtdirmsact,tbtdirmsdur]=crosscorel(filename,dataaligned,'all',0); 
% tbtcor=nanmedian(abs(tbtdircor))
% tbtcorstd=nanstd(abs(tbtdircor))
% corrcoef([tbtdirmsact tbtdirmsdur])

%% new code
%% behavioral analysis
% impact of previous trial's success/failure on current trial's
% success/failure and on saccade latency

behav=cellfun(@(x) x(1).trialnb, gsdata.allndata(~cellfun('isempty',gsdata.allsacdelay), 1),'UniformOutput',false); %NSS trials:
behav(:,2)=cellfun(@(x) [x(2:end).trialnb], gsdata.allndata(~cellfun('isempty',gsdata.allsacdelay), 1),'UniformOutput',false);%all SS trials
behav(:,3)=cellfun(@(x) x.nsst, gsdata.allsacdelay(~cellfun('isempty',gsdata.allsacdelay),1),'UniformOutput',false);%NSS RT
behav(:,4)=cellfun(@(x) x.ncst, gsdata.allsacdelay(~cellfun('isempty',gsdata.allsacdelay),1),'UniformOutput',false);%NCS RT

behav=behav(~cellfun('isempty',behav(:,2)),:);%removing bad apples

% Vincentizing not necessarily ideal
% see An evaluation of the Vincentizing method of forming group-level response time distributions.
% just pool RT distribution for respective conditions (trials preceded by
% either NSS or SS trial. Anderson-Darling test on distribution http://www.jaqm.ro/issues/volume-6,issue-3/pdfs/1_engmann_cousineau.pdf

behav(:,5)=cellfun(@(x,y,z) x(ismember(y,z+1)), behav(:,3),behav(:,1),behav(:,2),'UniformOutput',false);%RTs for NSS preceded by SS trial
behav(:,6)=cellfun(@(x,y,z) x(ismember(y,z+1)), behav(:,3),behav(:,1),behav(:,1),'UniformOutput',false);%RTs for NSS preceded by NSS trial

[H, adstat, critvalue] = adtest2([behav{:,5}], [behav{:,6}]); %no needt sort RTs

figure; subplot(1,2,1); hold on; hist([behav{:,6}]); hist([behav{:,5}]); 
h = findobj(gca,'Type','patch');
h(1).FaceColor = [0.8 0.2 0.1];
h(1).EdgeColor = 'w';
h(2).FaceColor = [0 0.5 0.5];
h(2).EdgeColor = 'w';
legend('RTs for NSS preceded by NSS trial','RTs for NSS preceded by SS trial')

[N,edges] = histcounts([behav{:,6}], 'Normalization', 'probability');
subplot(1,2,2); plot([edges(1)-diff(edges(1:2))/2 edges(2:end)-diff(edges)/2], cumsum([0 N]))
[N,edges] = histcounts([behav{:,5}], 'Normalization', 'probability');
hold on; plot([edges(1)-diff(edges(1:2))/2 edges(2:end)-diff(edges)/2], cumsum([0 N]))
axis('tight');box off;
title('Cumulative distribution')
legend('RTs for NSS preceded by NSS trial','RTs for NSS preceded by SS trial','location','southeast')


%% get indivudual trial's slope
% now with polynomial derivative, but hopefully with nonhomogeneous PP
% modeling in the future

% homogeneous poisson coeff
% y=poissrnd(0.5,1,100);
% y(y>0)=1;
% X=1:100; 
% coeff = glmfit(X, y,'poisson')
% 
% foo=gsdata.allgsndata{1, 1}(1).rast;
% foo=foo(:,200:600);
% [convsdf,convtrial]=conv_raster(foo,10);
% figure;plot(convsdf)
% 
% figure; hold on
% for tr=1:size(convtrial,1)
% %     plot(convtrial(43,:))
%     if mean(convtrial(tr,:))>0
%        bla(tr,1:2)= [tr  mean(convtrial(tr,:))];
%     end
% end
% 

% coeff = glmfit(1:81, foo(43,180:260),'poisson')
% yfit = glmval(coeff,1:81,-2);
% plot(1:81,yfit,'o')

p = polyfit(1:81, behav(43,180:260),4); %4 seems to work best for that sample 
k = polyder(p)    
plot(1:81,(1:81)*k(end-1)+k(end))

%get residuals as compared to gaussian smoothing 
bli=gauss_filtconv(behav(43,180:260),10); figure; hold on; plot(bli)