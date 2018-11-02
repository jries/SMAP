
function [locData,parameters,siteexplorer]=load_smlV3(filedat)
    locData=interfaces.LocalizationData;
    factorRefractiveMismatch=1.25;
    parameters=[];
    if isfield(filedat,'saveloc')
        saveloc=filedat.saveloc;

        locData.loc=saveloc.loc;
         if isfield(locData.loc,'z')
             locData.loc.znm=locData.loc.z*1000/factorRefractiveMismatch;
             disp('z corrected for refractive index mismatch and converted to nm')
         end

        s=length(saveloc.loc.xnm);
        if ~isfield(locData.loc,'locprecnm')
            if isfield(locData.loc,'locprecxnm')
                locData.loc.locprecnm=sqrt((locData.loc.locprecxnm.^2+locData.loc.locprecxnm.^2)/2);
            elseif isfield(locData.loc,'xerrnm')
                locData.loc.locprecnm=sqrt((locData.loc.xerrnm.^2+locData.loc.yerrnm.^2)/2);
            end
        end
        if ~isfield(locData.loc,'filenumber')
            locData.addloc('filenumber',ones(s(1),1,'uint16'));
        else
            locData.loc.filenumber=uint16(locData.loc.filenumber);
        end
        if ~isfield(saveloc.loc,'channel')
            locData.addloc('channel',zeros(s(1),1,'single'));
        end

        toremove={'xerrnm','yerrnm','bgerr','zerr'};
        locData.loc=myrmfield(locData.loc,toremove);

        if isfield(saveloc,'file')
            locData.files.file=saveloc.file;
            locData.files.filenumberEnd=length(saveloc.file);
            if length(locData.files.file)==1 %one file: overwrite file name
                locData.files.file.name=filedat.filename;
            end
        elseif isfield(saveloc,'info')      
            locData.files.filenumberEnd=1;
            locData.files.file.info=saveloc.info;
            if isfield(saveloc,'average')
                locData.files.file.average=saveloc.average;
            else
                locData.files.file.average=[];
            end
            locData.files.file.name=filedat.filename;
            locData.files.file.number=1;
            locData.files.file.numberOfTif=0;
            tif.image=[];
            tif.info=[];
            locData.files.file.tif=tif;
        end
        if isfield(saveloc,'siteexplorer')
            siteexplorer=saveloc.siteexplorer;
        else
            siteexplorer=[];
        end
        if isfield(saveloc,'history')
            locData.history=saveloc.history;
        end
        if isfield(saveloc,'fitparameters')
            locData.history{1}.children.fitparamters=saveloc.fitparameters;
        end
    end
    if isfield(filedat,'parameters')
        parameters=filedat.parameters;
    end

end
