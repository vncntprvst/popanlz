    
if strcmp(getenv('username'),'SommerVD') || strcmp(getenv('username'),'DangerZone')
    directory = 'C:\Data\Recordings\';
else
    directory = 'B:\data\Recordings\';
end
slash = '\';


% 
% algdir=[directory,'processed',slash,'aligned',slash];
% 
% load([algdir,,'_sac.mat']);     
dataaligned='';
%DPFl
 % 
% filename='R100L7A0_22800';
% filename='R103L7A1_23800';
% filename='S82L1A8_23301';

% filename='S99L2A5_10301';
% filename='R148L7A0_17503';
% filename='R142L7A1_16850';
% filename='S80L2A3_9901';

% filename='R116L6A1_14860';
filename='S79L3A3_9502';

[peakcct, peaksdf,tbtdircor,tbtdirmsact,tbtdirmsdur]=crosscorel(filename,dataaligned,'all',0); 
tbtcor=nanmedian(abs(tbtdircor))
tbtcorstd=nanstd(abs(tbtdircor))
corrcoef([tbtdirmsact tbtdirmsdur])
