global directory slash

if strcmp(getenv('username'),'SommerVD')
    directory = 'C:\Data\Recordings\';
elseif strcmp(getenv('username'),'DangerZone')
    directory = 'E:\data\Recordings\';
else
    directory = 'B:\data\Recordings\';
end
slash = '\';

%% get list of filename
postsacfiles={'R105L7A2_25160'
'R106L7A2_18600'
'R123L5P5_18400'
'R142L7A1_16850'
'R143L7A2_16900'
'R143L7A2_25651'
'R147L6P1_16200'
'R99L7A0_17400'
'S78L4A3_9651'
'S80L2A3_9901'
};
rampfiles={'R113L6A2_18900'
'R116L6A1_14860'
'R148L7A0_17503'
'R78L6A0_15240'
'S107L4A4_9151'
};
cat_files=postsacfiles;

%% get peak cross-correlation values: pretty reliable indicator to sort out presac, perisac and postsac activities
% possible limits are:
% pressac <-10ms before sac , >-10ms perisac <+10ms, <10ms postsac
peakcc=cell(length(cat_files),1);
for fnm=1:length(cat_files)
    peakcc{fnm}=crosscorel(cat_files{fnm},'active',0); %Get peakcc for all directions. Don't plot
end

%% area under curve: separate cells with low baseline FR and sharp burst from higher baseline neurons, especially ramping ones
% possible limit at 2000
auc=cell(length(cat_files),1);
for fnm=1:length(cat_files)
    auc{fnm}=findauc(cat_files{fnm},'all'); %Get auc for all directions
end
% get max auc
floor(cellfun(@(x) max(x), auc))