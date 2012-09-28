[sacnr,sacnrsessionlist,sacnrlocationlist,sacnrdepthlist,sacnrtasklist,compart,statscompile]=rec_class(2);

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

%all files from (crude classification) cortex, looking for strongest effect
cx_maxmdiff=maxmdiff(tcxidx | bcxidx);
cx_sacnr=sacnr(tcxidx | bcxidx);
[ord_cx_maxmdiff,smmdidx]=sort(cx_maxmdiff,'descend');
smmd_cx_sacnr=cx_sacnr(smmdidx);
%then top files - strongest effect
se_sacnr=smmd_cx_sacnr(1:5);

%best examples for each type of response 
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

