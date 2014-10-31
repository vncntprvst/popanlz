function [convsdf, convrasters, convrastsem]=conv_raster(rasters,fsigma,start,stop)
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

convrasters=NaN(size(rasters,1),stop-start+1);

for rast=1:size(rasters,1)
    convrasters(rast,:)=fullgauss_filtconv(rasters(rast,start:stop),fsigma,0).*1000;
end

convrasters=convrasters(:,fsigma*3+1:end-3*fsigma); %6 sigma ksize

convrastsem=std(convrasters)/ sqrt(size(convrasters,1));
convrastsem = convrastsem * 1.96; % 95% of the data will fall within 1.96 standard deviations of a normal distribution

convsdf=nanmean(convrasters);

% plot(convrasters,'k');

end              