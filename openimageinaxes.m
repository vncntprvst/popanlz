function openimageinaxes(imgname)

% open image
loadImage = importdata([imgname,'.png']); %dirfignames{1}
axes(findobj('Tag','figview')) % this gives focus to the axes
hold off; 
imshow(loadImage, [], 'InitialMagnification', 'fit');
end