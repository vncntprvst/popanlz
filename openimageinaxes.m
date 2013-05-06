function openimageinaxes(imgname,axeshandle)

% open image
loadImage = importdata([imgname{:},'.png']);
axes(axeshandle) 
set(axeshandle,'units','pixels');
figpos=get(axeshandle,'position');
figpos(3)=427;
figpos(4)=692; %figure size is 6920 x 4270 pixels
set(axeshandle,'position',figpos);
hold off; 
imshow(loadImage, [], 'InitialMagnification', 'fit');
axis(axeshandle, 'image')
end