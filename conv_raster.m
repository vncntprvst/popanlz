function convrasters=conv_raster(rasters,fsigma,start,stop)
%raster averaging

switch nargin
    case 1
        fsigma=10;
        start=1;
        stop=length(rasters);
    case 2
        start=1;
        stop=length(rasters);
end

convrasters=NaN(size(rasters,1),stop-start+2*fsigma+1);

for rast=1:size(rasters,1)
    convrasters(rast,:)=fullgauss_filtconv(rasters(rast,start-fsigma:stop+fsigma),fsigma,0);
end

convrasters=convrasters(:,fsigma+1:end-fsigma);
convrasters=nanmean(convrasters).*1000;
% plot(convrasters,'k');

end              