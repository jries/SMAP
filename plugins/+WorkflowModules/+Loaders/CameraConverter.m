classdef CameraConverter<interfaces.WorkflowModule
%     Interprets metadata and converts camera ADUs into photons- Metadata
%     can be overwritten manually or loaded from a SMAP _sml.mat data file.
    properties
        calfile='settings/CameraCalibration.xls';
        loc_cameraSettings=interfaces.metadataSMAP;
        loc_cameraSettingsStructure=struct('EMon',1,'emgain',1,'conversion',1,'offset',400,'cam_pixelsize_um',0.1,...
            'roi',[],'exposure',1,'timediff',0,'comment','','correctionfile','');        
        EMexcessNoise;
        calibrategain=false;
        calibratecounter
        calibrateimages
        offset
        adu2phot
        autocalhandle
        preview;
        gainmap;
        offsetmap;
       
    end
    methods
        function obj=CameraConverter(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.loc_cameraSettings=interfaces.metadataSMAP;
            fn=fieldnames(obj.loc_cameraSettingsStructure);
            for k=1:length(fn)
                obj.loc_cameraSettings.assigned.(fn{k})=false;
            end
            obj.propertiesToSave={'loc_cameraSettings'};
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function initGui(obj)
           initGui@interfaces.WorkflowModule(obj);
           obj.guihandles.loadmetadata.Callback={@loadmetadata_callback,obj};
           obj.guihandles.camparbutton.Callback={@camparbutton_callback,obj};
           obj.guihandles.calibrate.Callback={@calibrate_callback,obj};
            obj.outputParameters={'loc_cameraSettings'};
           obj.addSynchronization('loc_fileinfo_set',[],[],@obj.setmetadata)
           
        end
        function setmetadata(obj,overwrite)
            if nargin<2
                overwrite=false;
            end
            md=obj.loc_cameraSettings;
            
            lock=obj.getSingleGuiParameter('lockcampar');
            if lock
                disp('camera parameters locked, cannot be updated')
            end
            %update adu2phot
            if md.EMon
                emfactor=md.emgain;
            else
                emfactor=1;
            end
            md.pix2phot=md.conversion/emfactor;
            obj.adu2phot=md.pix2phot;
            obj.loc_cameraSettings=interfaces.metadataSMAP;
            obj.loc_cameraSettings=copyfields(obj.loc_cameraSettings,md);
            if ~overwrite 
                settings=obj.getPar('loc_fileinfo');
                fn=fieldnames(settings);
                if lock
                    lockedfields={'Width','Height','roi','exposure','emgain','EMon','conversion','offset','cam_pixelsize_um','timediff','comment','correctionfile'};
                    fn=setdiff(fn,lockedfields);
                end
                for k=1:length(fn)
                    if isfield(settings.assigned,fn{k}) && settings.assigned.(fn{k})
                        obj.loc_cameraSettings.(fn{k})=settings.(fn{k});
                    end
                end
            end
            
            obj.setPar('loc_cameraSettings',obj.loc_cameraSettings);
            obj.setPar('EMon',obj.loc_cameraSettings.EMon);
            obj.updatefileinfo;
%             if obj.loc_cameraSettings.EMon
%             	obj.EMexcessNoise=2;
%             else
%                 obj.EMexcessNoise=1;
%             end
        end
        function setcamerasettings(obj,fi)
            fn=fieldnames(fi);
            fn2=properties(obj.loc_cameraSettings);
            fna=intersect(fn,fn2);
            for k=1:length(fna)
                obj.loc_cameraSettings.(fna{k})=fi.(fna{k});
            end
        end
        function updatefileinfo(obj)
           fi=obj.getPar('loc_fileinfo');
           md2=obj.loc_cameraSettings;
           if isempty(fi)
               fn={};
           else
           fn=fieldnames(fi);
           end
           fn2=properties(md2);
           fna=intersect(fn,fn2);
           fi=copyfields(fi,md2,fna);
           obj.setPar('loc_fileinfo',fi);
        end
        function prerun(obj,p)
            pc=obj.loc_cameraSettings;
            if ~pc.EMon
                pc.emgain=1;
            end
            obj.offset=pc.offset;
            obj.adu2phot=(pc.conversion/pc.emgain);
            obj.preview=obj.getPar('loc_preview');
            loadcamcalibrationfile(obj);
        end
        function datao=run(obj,data,p)
           if isempty(data.data) %no image
               datao=data;
               return
           end
           imgp=makepositive(data.data);
           if p.correctcamera %apply offset and brightfield correction
               roi=p.loc_cameraSettings.roi;
               imphot=(single(imgp)-obj.offsetmap(roi(1)+1:roi(1)+roi(3),roi(2)+1:roi(2)+roi(4)))... %*obj.adu2phot...
                   .*obj.gainmap(roi(1)+1:roi(1)+roi(3),roi(2)+1:roi(2)+roi(4)); %XXX correct or exchanged???
               
           else
               imphot=(single(imgp)-obj.offset)*obj.adu2phot;
           end
               datao=data;
               datao.data=imphot;  
           if obj.preview && ~isempty(imphot)
               obj.setPar('preview_image',imphot);
           end
        end
       
    end
end


function out=makepositive(in)
if isa(in,'int16')
    out=single(in);
    out(out<0)=out(out<0)+2^16;
else
    out=in;
end
end

function calibrate_callback(a,b,obj)
cc=plugin('Analyze','calibrate','GainOffsetFFT');
cc.attachPar(obj.P);
if isempty(obj.autocalhandle)||~isvalid(obj.autocalhandle)
obj.autocalhandle=figure;
else
    delete(obj.autocalhandle.Children)
end
cc.handle=obj.autocalhandle;
cc.setGuiAppearence(struct('Vrim',100));
cc.makeGui;
f=obj.getPar('loc_filename');
if ~isempty(f)
    cc.loadfile(f);
% cc.setGuiParameters(struct('imagefile',f));
end
end

function camparbutton_callback(a,b,obj)
fn=fieldnames(obj.loc_cameraSettingsStructure);

fi=obj.getPar('loc_fileinfo');
obj.loc_cameraSettings=copyfields(obj.loc_cameraSettings,fi);
for k=length(fn):-1:1
    fields{k}=fn{k};
    if myisfield(fi,fn{k})
        defAns{k}=num2str(fi.(fn{k}));
    else
        defAns{k}=num2str(obj.loc_cameraSettings.(fn{k}));
    end
end
answer=inputdlg(fields,'Acquisition settings',1,defAns);
if ~isempty(answer) && ~ obj.getSingleGuiParameter('lockcampar')
    for k=1:length(fn)
        if isnumeric(obj.loc_cameraSettingsStructure.(fn{k}))||islogical(obj.loc_cameraSettingsStructure.(fn{k}))
            obj.loc_cameraSettings.(fn{k})=str2num(answer{k});
        else
            obj.loc_cameraSettings.(fn{k})=(answer{k});
        end
    end
obj.setmetadata(true);
end
if obj.getSingleGuiParameter('lockcampar')
    warning('cannot update camera paramters because they are locked')
end
end

function loadmetadata_callback(a,b,obj)

% return
finf=obj.getPar('loc_fileinfo');
if ~isempty(finf)
    ft=[(finf.basefile) filesep '*.*'];
else
    ft='*.*';
end

[f,p]=uigetfile(ft,'Select image (e.g. tiff) file or _sml.mat');
[~,~,ext ]=fileparts(f);
switch ext
    case '.mat'
        l=load([p f]);
        metadata=l.saveloc.file(1).info;
%     case '.txt'
%         par.metadatafile=[p f];
%         obj.setGuiParameters(par);
%         obj.updateGui;       
%         metadata=getmetadataMMtxt([p f]); 
    otherwise
        imloader=imageloaderAll([p f],finf,obj.P);
        metadata=imloader.getmetadata;
        metadata.allmetadata=metadata;
end
        obj.setcamerasettings(metadata);
        obj.updatefileinfo;
% obj.readmetadata;
end

function loadcamcalibrationfile(obj)
camfname=obj.loc_cameraSettings.correctionfile;
if strcmp(camfname(end-4:end),'.mat')
    camfname=camfname(1:end-5);
end
camfile=['settings' filesep 'cameras' filesep camfname '.mat'];
if exist(camfile,'file')
    l=load(camfile);
    
    [gainmap,offsetmap,varmap]=makegainoffsetCMOS(l,obj.loc_cameraSettings.exposure);
    
    obj.setPar('cam_varmap',single(varmap));
    
    obj.gainmap=gainmap;
    obj.offsetmap=offsetmap;
    
else
    if obj.getSingleGuiParameter('correctcamera')
        disp(['could not find camera correction file ' camfile])
    end
    obj.setPar('cam_varmap',[]);
end
end

function [gainmap,offsetmap,varmap]=makegainoffsetCMOS(l,exposuretime)
gainmap=l.gainmap;
offsetmap=l.offsetmap;
varmap=l.varmap;

varmap=varmap.*gainmap.^2; %???? in units of photons, to be used later directly
%XXXX here make gainmap, offsetmap
end

% function calculategain(img)
% img=double(img);
% m=mean(img,3);
% v=var(img,0,3);
% figure(88);
% plot(m(:),v(:),'.')
% gain=15.6;offs=200;em=200;
% pix2phot=gain/em;
% hold on
% plot(m(:),1/pix2phot*(m(:)-offs),'.')
% hold off
% end

function pard=guidef

pard.text.object=struct('Style','text','String','Metadata:');
pard.text.position=[1,1];
pard.text.Width=.7;
pard.text.Optional=true;
 
pard.loadmetadata.object=struct('Style','pushbutton','String','Load');
pard.loadmetadata.position=[1,1.6];
pard.loadmetadata.TooltipString=sprintf('Load camera settings metadata from image or _sml.mat file.');
pard.loadmetadata.Optional=true;
pard.loadmetadata.Width=0.4;
pard.calibrate.object=struct('Style','pushbutton','String','cal');
pard.calibrate.position=[1,2];
pard.calibrate.TooltipString=sprintf('calibrate gain and offset from images');
pard.calibrate.Optional=true;
pard.calibrate.Width=0.35;

pard.camparbutton.object=struct('Style','pushbutton','String','set Cam Parameters');
pard.camparbutton.position=[1,3.2];
pard.camparbutton.Width=1.3;
pard.camparbutton.TooltipString=sprintf('Edit camera acquisition parameters.');
pard.lockcampar.object=struct('Style','checkbox','String','Lock','Value',0);
pard.lockcampar.position=[1,4.5];
pard.lockcampar.Width=0.5;
pard.lockcampar.TooltipString=sprintf('Do not overwrite camera parameters automatically, but keep those set manually.');

pard.correctcamera.object=struct('Style','checkbox','String','Correct','Value',0);
pard.correctcamera.position=[1,2.6];
pard.correctcamera.Width=0.6;
pard.correctcamera.TooltipString=sprintf('Apply darkfield and brightfield correction.');


pard.plugininfo.type='WorkflowModule'; 
pard.plugininfo.description='Interprets metadata and converts camera ADUs into photons- Metadata can be overwritten manually or loaded from a SMAP _sml.mat data file.';
end