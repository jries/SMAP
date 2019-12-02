function filestruct=locsaveFileinfo(obj)
filestruct=initfile;
% tifs=struct('image',[],'info',[]);
finfo=prop2struct(obj.getPar('loc_fileinfo'));
sfile=finfo.basefile;

%single image files

if isempty(sfile)
    sfile=filestruct.name;
end


if isempty(strfind(sfile,'_sml.mat'))
    [path,filen]=fileparts(sfile);
    filename=[path filesep filen '_sml.mat']; 
else
    filename=sfile;
end
infost=prop2struct(obj.getPar('loc_cameraSettings'));
infost=copyfields(infost,finfo);

newpix=obj.getPar('overwrite_pixelsize');
if ~isempty(newpix) 
    infost.cam_pixelsize_um=newpix;
end
filestruct.info=infost;
filestruct.name=filename;
filestruct.number=obj.filenumber;
% filestruct=struct('info',infost,'average',[],'name',filename,...
%     'number',obj.filenumber,...
%     'numberOfTif',0,'tif',tifs,'raw',struct('image',[],'frame',[]));
if ~myisfield(filestruct.info,'roi')||isempty(filestruct.info.roi)
    try
        mm=imfinfo(finfo.filename);
        obj.fileinfo.roi=[0 0 mm(1).Width mm(1).Height];

    catch
        try
            filestruct.info.roi=[0 0 filestruct.info.Width filestruct.info.Height];
        catch
            filestruct.info.roi=[0 0 512 512];
        end
    end
    filestruct.info.roi=obj.fileinfo.roi;
end
% globT=obj.getPar('loc_globaltransform');
% if isa(globT,'interfaces.LocTransformN') && ~isempty(globT) && myisfield(globT,'info')
%     tinfoh=globT.info{1};
%     filestruct.info.roi(3)=min(filestruct.info.roi(3),tinfoh.xrange(2));filestruct.info.roi(4)=min(filestruct.info.roi(4),tinfoh.yrange(2));
% end
