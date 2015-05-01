function [convsdf, convrasters, convrastsem]=conv_raster(rasters,fsigma,start,stop,meanFR,stdFR)
%raster averaging

switch nargin
    case 1
        fsigma=10;
        start=1;
        stop=size(rasters,2);
        meanFR=0;
        stdFR=1;
    case 2
        start=1;
        stop=size(rasters,2);
        meanFR=0;
        stdFR=1;
    case 4
        meanFR=0;
        stdFR=1;
end

convrasters=NaN(size(rasters,1),stop-start+1);

for rast=1:size(rasters,1)
    convrasters(rast,:)=fullgauss_filtconv(rasters(rast,start:stop),fsigma,0).*1000;
%     figure; plot(convrasters(rast,:));
end

convrasters=convrasters(:,fsigma*3+1:end-3*fsigma); %6 sigma ksize

%normalization
convrasters=(convrasters-meanFR)./stdFR;

%SEM
convrastsem=std(convrasters)/ sqrt(size(convrasters,1));
convrastsem = convrastsem * 1.96; % 95% of the data will fall within 1.96 standard deviations of a normal distribution

convsdf=nanmean(convrasters);

% figure; plot(convrasters,'k');

end              