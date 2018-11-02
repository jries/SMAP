function [locData,parameters,siteexplorer]=loadfitposV2(filedat)
locData=interfaces.LocalizationData;
factorRefractiveMismatch=1.25
average=[];
%     load(file)
fitpos=filedat.fitpos;
infox=filedat.infox;

    s=size(fitpos);
%     locData.info=infox;
    locData.addloc('frame',fitpos(:,1));
    pixsizenm=infox.pixsize*1000;
    locData.setloc('xnm',fitpos(:,3)*pixsizenm);
    locData.setloc('ynm',fitpos(:,2)*pixsizenm);
    locData.setloc('phot',fitpos(:,4));
    locData.setloc('bg',fitpos(:,5));
    if min(fitpos(:,6))<0 %z-loc, test better in the future
        locData.setloc('znm',fitpos(:,6)*1000/factorRefractiveMismatch);
        locData.setloc('PSFxnm',fitpos(:,6)*0+150);
        locData.setloc('locprecznm',fitpos(:,13)*1000/factorRefractiveMismatch)
    else
        locData.setloc('PSFxnm',fitpos(:,6)*pixsizenm);
        locData.setloc('PSFynm',fitpos(:,7)*pixsizenm);
    end
    locData.setloc('locprecnm',fitpos(:,19)*pixsizenm);
    if s(2)>24
    locData.setloc('channel',fitpos(:,25));
%     locData.addloc('original_channel',fitpos(:,25));
    else
        locData.setloc('channel',zeros(s(1),1,'single'));
    end
    locData.setloc('filenumber',zeros(s(1),1,'uint8')+1);
    locData.setloc('logLikelihood',fitpos(:,16));
    
    if s(2)>30
        locData.setloc('int1',fitpos(:,31));
        locData.setloc('int2',fitpos(:,33));
    end
    
    infox.roi=[0 0 infox.Width infox.Height];
    
    locData.files.filenumberEnd=1;
    locData.files.file.info=infox;
        locData.files.file.average=average;
    locData.files.file.name=filedat.filename;
    locData.files.file.number=1;
    locData.files.file.numberOfTif=0;
    tif.image=[];
    tif.info=[];
    locData.files.file.tif=tif;
    
    %sometimes ROI is wrong
    parameters=[];
    siteexplorer=[];
end