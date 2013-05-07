
subject={'Rigel','Sixx','Hilda'};
%emdirections={'upward','up_left','leftward','down_left','downward','down_right','rightward','up_right','all'};
emdirections={'all'};
for monknum=1:3
% alldata=findssrt_batch(monknum,emdirections);
alldata=SSRT_TachoMP(monknum,emdirections);
meansaclat=nan(length(emdirections),1);
for andir=1:length(emdirections)
dalldata=alldata(arrayfun(@(x) strcmp(x.dir,emdirections{andir}), alldata));
dnssdelay=cellfun(@(x) x.all, {dalldata(~cellfun('isempty',{dalldata.NSSsacdelay})).NSSsacdelay}, 'UniformOutput', false);
meansaclat(andir)=round(mean([dnssdelay{:}]));
meansaclatdisp=sprintf('NSS mean sac latency in %s direction is: %d', emdirections{andir},meansaclat(andir));
disp(meansaclatdisp);
end

% display saccade latency 
polx=0:2*pi/8:2*pi; %nine values to close the loop
meansaclat=meansaclat([7 8 1:6 7]); %redistribute from rightward to rightward, counterclockwise

% Create polar plot
figure('Color',[1 1 1]);
polar(polx,meansaclat');
saclatlines=findobj(gcf,'type','line');
set(saclatlines(end-1),'LineWidth',3,'color',[0.2 0.5 0.8]);
title({'Mean saccade latency (ms) for no-stop signal trials, by direction'},...
    'FontSize',14,...
    'FontName','Calibri');
legend1 = legend(gca,subject{monknum});
set(legend1,...
    'Position',[0.829671448813315 0.590575396825395 0.0892667375132837 0.0505952380952381],...
    'FontSize',12,...
    'FontName','Calibri');

end
