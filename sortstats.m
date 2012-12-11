if strcmp(getenv('username'),'SommerVD')
    directory = 'C:\Data\Recordings\';
elseif strcmp(getenv('username'),'DangerZone')
    directory = 'E:\data\Recordings\';
else
    directory = 'B:\data\Recordings\';
end
slash = '\';

%% uniquefiles
uniquefiles={'R100L7A0_2280'
'R100L7A0_2300'
'R100L7A0_2410'
'R102L7A1_1750'
'R102L7A1_2420'
'R103L7A1_2380'
'R104L7A2_2020'
'R104L7A2_2080'
'R105L7A2_2516'
'R106L7A2_1765'
'R106L7A2_1810'
'R106L7A2_1860'
'R107L6A3_1851'
'R107L6A3_1917'
'R107L6A3_2097'
'R109L6A1_1640'
'R109L6A1_1700'
'R110L6A2_1875'
'R110L6A2_1980'
'R110L6A2_1986'
'R110L6A2_2010'
'R111L6A2_2311'
'R112L6A2_1961'
'R112L6A2_1985'
'R112L6A2_2045'
'R112L6A2_2066'
'R113L6A2_1890'
'R113L6A2_2231'
'R113L6A2_2525'
'R144L7A2_1900'
'R25L5A3_1400'
'R29L7A1_1750'
'R32L6P1_1050'
'R34L6A0_2000'
'R36L6A0_1560'
'R37L6A0_2150'
'R39L7A0_1800'
'R42L7A1_1850'
'R43L7A1_1850'
'R44L7A1_2170'
'R61L6A0_1492'
'R61L6A0_1880'
'R70L6A0_1600'
'R70L6A0_1826'
'R70L6A0_1840'
'R70L6A0_1896'
'R71L6A0_1649'
'R71L6A0_1673'
'R72L6A0_1558'
'R72L6A0_1580'
'R74L6A0_1662'
'R74L6A0_1690'
'R75L6A0_1658'
'R93L6A1_2220'
'R94L6A3_1810'
'R94L6A3_2100'
'R95L6A2_1920'
'R95L6A2_1940'
'R96L7A1_1910'
'R96L7A1_2000'
'R97L7A1_1790'
'R98L7A0_1860'
'R99L7A0_1880'}

%% call recordings classification
%[sacnr,sacnrsessionlist,sacnrlocationlist,sacnrdepthlist,sacnrtasklist,compart,statscompile]=rec_class(2);
[rsacnr,rsacnrsessionlist,rsacnrlocationlist,rsacnrdepthlist,rsacnrtasklist,rcompart,rstatscompile]=rec_class(1);
[ssacnr,ssacnrsessionlist,ssacnrlocationlist,ssacnrdepthlist,ssacnrtasklist,scompart,sstatscompile]=rec_class(2);
sacnr=[rsacnr;ssacnr];
sacnrsessionlist=[rsacnrsessionlist;ssacnrsessionlist];
sacnrlocationlist=[rsacnrlocationlist;ssacnrlocationlist];
sacnrdepthlist=[rsacnrdepthlist;ssacnrdepthlist];
sacnrtasklist=[rsacnrtasklist;ssacnrtasklist];
compart=[rcompart;scompart];
% add empty fields to sstatcompile to concatenate with rstatcompile
for i=1:length(sstatscompile)
    emptystruct(i,1:size(rstatscompile,2)-1)=struct('p',[],'h',[],'bestmeandiff',[],'recrasters',[]);
end
sstatscompile=[sstatscompile emptystruct];
statscompile=[rstatscompile;sstatscompile];

rigelfiles=~cellfun(@isempty,regexp(sacnr,'^R'));
sixxfiles=~cellfun(@isempty,regexp(sacnr,'^S'));
tcxidx=~cellfun(@isempty,regexp(compart,'top_cortex'));
cdnidx=~cellfun(@isempty,regexp(compart,'dentate'));
bcxidx=~cellfun(@isempty,regexp(compart,'bottom_cortex'));

maxmdiff=nan(size(statscompile));
effectype=cell(size(statscompile));
effectcategories={'sactobase','pre_post','sharp'};
for i=1:size(statscompile,1)
    maxmdiff(i)=max(cellfun(@(x) round(max(x)), statscompile(i,1).bestmeandiff));
    treffectype=[statscompile(i,1).h{:}];
    if length(treffectype)>1
        treffectype=[treffectype{:}];
        summeffect=[0,0,0];
        for j=1:3:(length(treffectype))
            summeffect=summeffect+treffectype(j:j+2);
        end
        effectype(i,1:sum(logical(summeffect)))=effectcategories(logical(summeffect));
    else
        effectype(i,1:sum(logical(treffectype{:})))=effectcategories(logical([treffectype{:}]));
    end
end
%% writing rough loc and max mean diff to xls file
%packaging data
rmaxmdiff=maxmdiff(rigelfiles);
smaxmdiff=maxmdiff(sixxfiles);

% get number of row in "database"
exl = actxserver('excel.application');
exlWkbk = exl.Workbooks;
exlFile = exlWkbk.Open([directory 'procdata.xlsx']);
for monknum=1:2
    exlSheet(monknum) = exlFile.Sheets.Item(monknum);% e.g.: 2 = Sixx
    robj(monknum) = exlSheet(monknum).Columns.End(4);
    numrows(monknum) = robj(monknum).row;
    if numrows(monknum)==1048576 %empty document
        numrows(monknum)=1;
    end 
end
Quit(exl);

for monknum=1:2
    [~,pfilelist] = xlsread([directory,'procdata.xlsx'],monknum,['A2:A' num2str(numrows(monknum))]);
    if monknum==1
        for mkfl=1:size(rsacnr,1)
            wline=find(ismember(pfilelist,rsacnr(mkfl)))+1;
            xlswrite([directory,'procdata.xlsx'], rcompart(mkfl), monknum, sprintf('H%d',wline));
            xlswrite([directory,'procdata.xlsx'], {rmaxmdiff(mkfl)}, monknum, sprintf('M%d',wline));
        end
    elseif monknum==2
        for mkfl=1:size(ssacnr,1)
            wline=find(ismember(pfilelist,ssacnr(mkfl)))+1;
            xlswrite([directory,'procdata.xlsx'], scompart(mkfl), monknum, sprintf('H%d',wline));
            xlswrite([directory,'procdata.xlsx'], {smaxmdiff(mkfl)}, monknum, sprintf('M%d',wline));
        end
    elseif monknum==3
        for mkfl=1:size(hsacnr,1)
            wline=find(ismember(pfilelist,hsacnr(mkfl)))+1;
            xlswrite([directory,'procdata.xlsx'], {hcompart(mkfl)}, monknum, sprintf('H%d',wline));
            xlswrite([directory,'procdata.xlsx'], {hmaxmdiff(mkfl)}, monknum, sprintf('M%d',wline));
        end
    end
end


%% continue analysis on cortex files
%all files from (crude classification) cortex, looking for strongest effect
cx_maxmdiff=maxmdiff(tcxidx | bcxidx);
cx_sacnr=sacnr(tcxidx | bcxidx);
[smmd_cx_maxmdiff,smmdidx]=sort(cx_maxmdiff,'descend');
smmd_cx_sacnr=cx_sacnr(smmdidx);

%applying cortex treatment to tasklist
cx_sacnrtasklist=sacnrtasklist(tcxidx | bcxidx);
smmd_cx_sacnrtasklist=cx_sacnrtasklist(smmdidx);

%'top cx' data
tcx_sacnr=sacnr(tcxidx);
tcx_maxmdiff=maxmdiff(tcxidx);
[smmd_tcx_maxmdiff,smmdidx]=sort(tcx_maxmdiff,'descend');
smmd_tcx_sacnr=tcx_sacnr(smmdidx);
tcxdata=[smmd_tcx_sacnr num2cell(smmd_tcx_maxmdiff)];

%'bottom cx' data
bcx_sacnr=sacnr(bcxidx);
bcx_maxmdiff=maxmdiff(bcxidx);
[smmd_bcx_maxmdiff,smmdidx]=sort(bcx_maxmdiff,'descend');
smmd_bcx_sacnr=bcx_sacnr(smmdidx);
bcxdata=[smmd_bcx_sacnr num2cell(smmd_bcx_maxmdiff)];

%top files - strongest effect
se_sacnr=smmd_cx_sacnr(1:5);

%best examples for each type of response (effectype calculated based
%bestmeandiff)
stbeff=sum(cellfun(@(x) strcmp([x],'sactobase'), effectype),2);
ppeff=sum(cellfun(@(x) strcmp([x],'pre_post'), effectype),2);
sharpeff=sum(cellfun(@(x) strcmp([x],'sharp'), effectype),2);
%narrowing down to cortex files
cx_stbeff=stbeff(tcxidx | bcxidx);
cx_ppeff=ppeff(tcxidx | bcxidx);
cx_sharpeff=sharpeff(tcxidx | bcxidx);

onlystb=cx_stbeff & (~(cx_ppeff | cx_sharpeff));
onlypp= cx_ppeff & (~(cx_stbeff | cx_sharpeff));
onlysharp=cx_sharpeff & (~(cx_ppeff | cx_stbeff));

stb_pp= (cx_stbeff & cx_ppeff) & (~(cx_sharpeff));
stb_sharp= (cx_stbeff & cx_sharpeff) & (~(cx_ppeff));
pp_sharp= (cx_ppeff & cx_sharpeff) & (~(cx_stbeff));

alleffects= cx_stbeff & cx_ppeff & cx_sharpeff;
%finally find best files for each effect category
smmd_onlystb=smmd_cx_sacnr(onlystb(smmdidx));
smmd_onlypp=smmd_cx_sacnr(onlypp(smmdidx));
smmd_onlysharp=smmd_cx_sacnr(onlysharp(smmdidx));
smmd_stb_pp=smmd_cx_sacnr(stb_pp(smmdidx));
smmd_stb_sharp=smmd_cx_sacnr(stb_sharp(smmdidx));
smmd_pp_sharp=smmd_cx_sacnr(pp_sharp(smmdidx));
smmd_alleffects=smmd_cx_sacnr(alleffects(smmdidx));



%% plotting files by category
% only sac to baseline
effectcat='onlystb';
SummaryPlot(effectcat,smmd_onlystb,smmd_cx_sacnrtasklist(onlystb(smmdidx)));
close all force;
% only pre post sac
effectcat='prepost';
SummaryPlot(effectcat,smmd_onlypp,smmd_cx_sacnrtasklist(onlypp(smmdidx)));
close all force;
% only sharp sac
effectcat='sharpsac';
SummaryPlot(effectcat,smmd_onlysharp,smmd_cx_sacnrtasklist(onlysharp(smmdidx)));
close all force;
% sac to baseline and prepost
effectcat='stbpp';
SummaryPlot(effectcat,smmd_stb_pp,smmd_cx_sacnrtasklist(stb_pp(smmdidx)));
close all force;
% sac to baseline and sharp
effectcat='stbsh';
SummaryPlot(effectcat,smmd_stb_sharp,smmd_cx_sacnrtasklist(stb_sharp(smmdidx)));
close all force;
% prepost and sharp
effectcat='ppsh';
SummaryPlot(effectcat,smmd_pp_sharp,smmd_cx_sacnrtasklist(pp_sharp(smmdidx)));
close all force;
% all effects
effectcat='alleffects';
SummaryPlot(effectcat,smmd_alleffects,smmd_cx_sacnrtasklist(alleffects(smmdidx)));
close all force;

% over 19 max mean diff
effectcat='20over';
SummaryPlot(effectcat,smmd_cx_sacnr(find(smmd_cx_maxmdiff>20)),smmd_cx_sacnrtasklist(find(smmd_cx_maxmdiff>20)));
close all force;

%just looking at stats from statscompile:
% filetoget=smmd_onlystb{1};
% ismember(filetoget,sacnr);
% [~,fidx]=ismember(filetoget,sacnr);
% statscompile(fidx).p{:}{:};
% statscompile(fidx).h{:}{:};



