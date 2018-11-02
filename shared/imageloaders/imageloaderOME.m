classdef imageloaderOME<interfaces.imageloaderSMAP
    %imageloaderMM image loader for micromanager single tiff files
    %   Detailed explanation goes here
    
    properties
        calfile='settings/CameraCalibration.xls';
        reader
        seriesnumber=0;
%         allmetadatatags
    end
    
    methods
        function obj=imageloaderOME(varargin)
            obj@interfaces.imageloaderSMAP(varargin{:});
        end
        function openi(obj,file)
            obj.file=file;
            obj.reader=bfGetReader(file);
            md=obj.getmetadata;
            
        end
        function closei(obj)
            obj.reader.close;
        end
%         function mdo=getmetadata(obj)
% %             mdo=getmetadataome(obj);
%             mdo=getmetadatacam(obj);
%              
%         end
        function image=getimagei(obj,frame)
            image=readstack(obj,frame);
        end
        
        function allmd=getmetadatatagsi(obj)
        [ph,fh,ext]=fileparts(obj.file);
        exclude={};
        switch ext
            case '.nd2' %Nikon
                sm=obj.reader.getGlobalMetadata;
                exclude={'X position for position','Y position for position','Z position for position','PFS Status #','PFS Offset #'};
                parselist=true;
            case '.lif'
                setseries(obj);
                cm=obj.reader.getCoreMetadataList;
                cm1=cm.get(0);
                sm=cm1.seriesMetadata;
                parselist=true;
            case '.tif'
                parselist=true;
                   sm=obj.reader.getGlobalMetadata;
            otherwise 
                parselist=true;
                sm=obj.reader.getGlobalMetadata;
                
        end
                
        if parselist  &&      ~isempty(sm) 
                k=sm.keys;
                ind=1;
                allmd=[];
                while (ismethod(k,'hasNext') && k.hasNext) || (ismethod(k,'hasMoreElements') && k.hasMoreElements) 
                    kh=k.nextElement;
                    if ~any(strncmp(exclude,kh,10))
                        try
                            v=sm.get(kh);
                            if ~isempty(v)
                                allmd(ind,1:2)={kh ,v};
                                ind=ind+1;
                            end
                        catch
                        end
                     
                    end
                end
                allome=getmetadatatagsome(obj);
                allmd=vertcat(allmd,allome);
        else
            allmd=[];
        end
        try
        allcore=getallcoremtadata(obj.reader.getCoreMetadataList);
        allmd=vertcat(allmd,allcore);

        catch
        end
        f=getnumberofframes(allmd);
        allmd(end+1,:)={'frames direct',num2str(f)};        
        obj.allmetadatatags=allmd;
        end
        
    end
    
end

function f=getnumberofframes(allmd)
numf=contains(allmd(:,1),{'sizeZ','sizeC','sizeT','imageCount','Positions'});
f=nanmax(str2double(allmd(numf,2)));
% obj.metadata.allmetadata(end+1,:)={'frames direct',num2str(f)};
end

function out=getallcoremtadata(cm)
sm=cm.get(0);
list={'sizeX','sizeY','sizeZ','sizeC','sizeT','pixelType','bitsPerPixel','imageCount','dimensionOrder'};
for k=1:length(list)
    try
        out{k,1}=list{k};
        v=(sm.(list{k}));
        if isnumeric(v)
            v=num2str(v);
        else
            v=char(v);
        end
        out{k,2}=v;
    end
end
end
% out(1,1:2)={'sizeX',sm.sizeX
% end
function metao=getmetadataome(obj)
reader=obj.reader;
omemeta=reader.getMetadataStore;

[ph,fh,ext]=fileparts(obj.file);
fn={};
indseries=0;
switch ext
    case '.nd2' %Nikon
        meta=getMetaNd2(reader);
        fn=fieldnames(meta);
    case '.lif'
        meta=getMetaLif(reader);
        fn=fieldnames(meta);
        %determine series
        seri=reader.getSeriesCount;
        ind=1;
        series=[];
        message={};
        for k=seri:-1:1
            reader.setSeries(k-1);
            numim(k)=reader.getImageCount;
            if numim(k)>1
                series(ind)=k;
                message{ind}=['S' num2str(k) ', ' num2str(numim(k)) ' frames'];
                ind=ind+1;
            end
        end
        if ind>2
            if isempty(obj.seriesnumber)
                selected=listdlg('ListString',message,'SelectionMode','single','Name','Select data set');
                largeseries=series(selected);
                 obj.seriesnumber=largeseries-1;
            else
                largeseries=obj.seriesnumber+1;
            end
        else
            [~,largeseries]=max(numim);
        end
        reader.setSeries(largeseries-1);
        indseries=largeseries-1;
    otherwise
        meta=[];
end
try
m2.cam_pixelsize_um=[double(omemeta.getPixelsPhysicalSizeX(indseries).value()) double(omemeta.getPixelsPhysicalSizeY(indseries).value())];
fn=[fn fieldnames(m2)];
catch
    m2=[];
end
obj.metadata=copyfields(obj.metadata,meta);
obj.metadata=copyfields(obj.metadata,m2);
obj.metadata.Width=double(omemeta.getPixelsSizeX(indseries).getValue());
obj.metadata.Height=double(omemeta.getPixelsSizeY(indseries).getValue());
obj.metadata.numberOfFrames=max(double(omemeta.getPixelsSizeT(indseries).getValue()),double(omemeta.getPixelsSizeZ(indseries).getValue()));
obj.metadata.basefile=[ph filesep fh];
obj.metadata.roi=[0 0 obj.metadata.Width obj.metadata.Height];

fn=[makehorz(fn) makehorz({'Width','Height','numberOfFrames','basefile'})];

for k=1:length(fn)
    obj.metadata.assigned.(fn{k})=true;
end
metao=obj.metadata;
end

function setseries(obj)
reader=obj.reader;
seri=reader.getSeriesCount;
        ind=1;
        series=[];
        message={};
        for k=seri:-1:1
            reader.setSeries(k-1);
            numim(k)=reader.getImageCount;
            if numim(k)>1
                series(ind)=k;
                message{ind}=['S' num2str(k) ', ' num2str(numim(k)) ' frames'];
                ind=ind+1;
            end
        end
        if ind>2
            if isempty(obj.seriesnumber)
                selected=listdlg('ListString',message,'SelectionMode','single','Name','Select data set');
                largeseries=series(selected);
                 obj.seriesnumber=largeseries-1;
            else
                largeseries=obj.seriesnumber+1;
            end
        else
            [~,largeseries]=max(numim);
        end
        reader.setSeries(largeseries-1);
        obj.seriesnumber=largeseries-1;
end



function md=getMetaLif(reader)

cm=reader.getCoreMetadataList;
cm1=cm.get(1);
sm=cm1.seriesMetadata;
% md.emgain=str2double(sm.get('ProcessingHistory|ATLCameraSettingDefinition|EMGainValue'));
% md.conversion=str2double(sm.get('ProcessingHistory|ATLCameraSettingDefinition|GainValue'));
% md.EMon=str2double(sm.get('ProcessingHistory|ATLCameraSettingDefinition|CanDoEMGain'));
% md.exposure=1000*str2double(sm.get('ProcessingHistory|ATLCameraSettingDefinition|WideFieldChannelConfigurator|SameExposureTime'));
md.emgain=str2double(sm.get('ATLCameraSettingDefinition|EMGainValue'));
md.conversion=str2double(sm.get('ATLCameraSettingDefinition|GainValue'));
md.EMon=str2double(sm.get('ATLCameraSettingDefinition|CanDoEMGain'));
md.exposure=1000*str2double(sm.get('ATLCameraSettingDefinition|WideFieldChannelConfigurator|SameExposureTime'));
md.timediff=md.exposure;

fn=fieldnames(md);
for k=1:length(fn)
    if isnan(md.(fn{k}))
        md=rmfield(md,fn{k});
    end
end

k=sm.keys;
ind=1;
while k.hasNext
    kh=k.nextElement;
    allmd{ind}=([kh ',' sm.get(kh)]);
    ind=ind+1;

end
md.allmetadata.omeLif=allmd;
end



function allmd=getmetadatatagsome(obj)
allmd=[];
try
    omemeta=obj.reader.getMetadataStore;
    allmd(1,:)={'getPixelsPhysicalSizeX',double(omemeta.getPixelsPhysicalSizeX(obj.seriesnumber).value())};
    allmd(2,:)={'getPixelsSizeX',double(omemeta.getPixelsSizeX(obj.seriesnumber).getValue())};
    allmd(3,:)={'getPixelsSizeY',double(omemeta.getPixelsSizeY(obj.seriesnumber).getValue())};
    allmd(4,:)={'getPixelsSizeT',double(omemeta.getPixelsSizeT(obj.seriesnumber).getValue())};
    allmd(5,:)={'getPixelsSizeZ',double(omemeta.getPixelsSizeZ(obj.seriesnumber).getValue())};  
catch
    disp('could not get core metadata')
end
end

function meta=getMetaNd2(reader)
    m=reader.getGlobalMetadata;
    meta.EMon=str2num(m.get('EnableGainMultiplier'));
    meta.exposure= str2num(m.get('Exposure'));
    meta.emgain=str2num(m.get('GainMultiplier'));
    meta.conversion=str2num(m.get('ConversionGain'));
    meta.timediff=meta.exposure;
end




function image=readstack(obj,imagenumber)
   if imagenumber<=obj.reader.getImageCount()
       image=bfGetPlane(obj.reader,imagenumber);
   else
       image=[];
   end
end

