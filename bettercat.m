global directory slash

if strcmp(getenv('username'),'SommerVD') || strcmp(getenv('username'),'DangerZone')
    directory = 'C:\Data\Recordings\';
else
    directory = 'B:\data\Recordings\';
end
slash = '\';

%% get peak cross-correlation values
cat_files={'R99L7A0_17400'
'S78L4A3_9651'
'S80L2A3_9901'
};
%ramp
% {'R113L6A2_18900'
% 'R116L6A1_14860'
% 'R148L7A0_17503'
% 'R78L6A0_15240'
% 'S107L4A4_9151'
% };
peakcc=cell(length(cat_files),1);
for fnm=1:length(cat_files)
    peakcc{fnm}=crosscorel(cat_files{fnm},0); %don't plot
end