[sacnr,sacnrsessionlist,sacnrlocationlist,sacnrdepthlist,sacnrtasklist,compart,statscompile]=rec_class(2);
maxmdiff=nan(size(statscompile));
for i=1:size(statscompile,1)
    maxmdiff(i)=max(cellfun(@(x) round(max(x)), statscompile(i,1).bestmeandiff));
end

tcxidx=~cellfun(@isempty,regexp(compart,'top_cortex'));
cdnidx=~cellfun(@isempty,regexp(compart,'dentate'));
bcxidx=~cellfun(@isempty,regexp(compart,'bottom_cortex'));

cx_maxmdiff=maxmdiff(tcxidx | bcxidx);
cx_sacnr=sacnr(tcxidx | bcxidx);
[ord_cx_maxmdiff,smmdidx]=sort(cx_maxmdiff,'descend');
smmd_cx_sacnr=cx_sacnr(smmdidx);

%top file - strongest effect
se_sacnr=smmd_cx_sacnr(1:5);

%best examples for each type of response 