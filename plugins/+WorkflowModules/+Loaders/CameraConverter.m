classdef CameraConverter<interfaces.WorkflowModule
%     Interprets metadata and converts camera ADUs into photons- Metadata
%     can be overwritten manually or loaded from a SMAP _sml.mat data file.
    properties
        calfile;
        loc_cameraSettings=interfaces.metadataSMAP;
        loc_cameraSettingsStructure=struct('EMon',1,'emgain',1,'conversion',1,'offset',400,'cam_pixelsize_um',0.1,...
            'roi',[],'exposure',1,'timediff',0,'comment','','correctionfile','','imagemetadata','');        
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
        varmap
        offsetmapuse;
        gainuse;
        scmosroi;
        rawaverage
        rawaveragecounter
        rawaveragetype
        rawimagestruct
        rawimagecounter
       
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
            pard=guidef(obj);
        end
        function initGui(obj)
           initGui@interfaces.WorkflowModule(obj);
           obj.guihandles.loadmetadata.Callback={@loadmetadata_callback,obj};
           obj.guihandles.camparbutton.Callback={@camparbutton_callback,obj};
           obj.guihandles.calibrate.Callback={@calibrate_callback,obj};
            obj.outputParameters={'loc_cameraSettings'};
           obj.addSynchronization('loc_fileinfo_set',[],[],@obj.setmetadata)
            obj.inputParameters={'diffrawframes'};          
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
                    lockedfields={'Width','Height','roi','exposure','emgain','EMon','conversion','offset','cam_pixelsize_um','timediff','comment','correctionfile','imagemetadata'};
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
            mirrorem_callback(0,0,obj)
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
            
            obj.rawimagecounter=1;
            obj.rawaveragecounter=1;
            obj.rawaverage=[];
            obj.offsetmapuse=[];
            mirrorem_callback(0,0,obj);
            obj.setPar('cam_varmap',[]);
               
            
            %             if fileinf.EMon && p.mirrorem  %if em gain on and mirrorem on: switch roi
%                 %It seems that on the Andor the roi is independent on the
%                 %mode, 
% %                 if any(fileinf.roi(1:2)>0) %if roi(1:2)=[0 0] it is likely that roi was not read out and set to default.
% %                     fileinf.roi(1)=512-fileinf.roi(1)-fileinf.roi(3);
% %                 end
%                 fileinf.EMmirror=true;
%             else 
%                 fileinf.EMmirror=false;
%             end
%             obj.mirrorem=fileinf.EMmirror;
        end
        function datao=run(obj,data,p)
            if data.eof %transmit image stack
                obj.rawimagestruct(1).image=cast(obj.rawaverage/(obj.rawaveragecounter-1),'like', obj.rawaveragetype);
                obj.rawimagestruct(1).frame=0;
                
                obj.setPar('rawimagestack',obj.rawimagestruct(1:obj.rawimagecounter));
            end
           if isempty(data.data) %no image
               datao=data;
               return
           end
           
           img=data.data;
           if p.correctcamera %apply offset and brightfield correction
               if isempty(obj.offsetmapuse)
                    loadcamcalibrationfile(obj,p,img);
               end
           end

%            if iscell(img)
%                for k=length(img):-1:1
%                    imphot{k}=convertsingleimage(obj,img{k},p,obj.offsetmapuse);
%                end
%            else
               imphot=convertsingleimage(obj,img,p,obj.offsetmapuse);
%            end
              
%                imphot=(single(imgp)-obj.offsetmapuse)*obj.gainuse;   
% %                imphot=(single(imgp)-obj.offsetmap(roi(1)+1:roi(1)+roi(3),roi(2)+1:roi(2)+roi(4)))... %*obj.adu2phot...
% %                    ./obj.gainmap(roi(1)+1:roi(1)+roi(3),roi(2)+1:roi(2)+roi(4)); %XXX correct or exchanged??? a: looks fine
%            else
%            imgp=makepositive(data.data);
%            if p.emmirror && obj.loc_cameraSettings.EMon  
%                     imgp=imgp(:,end:-1:1);
%            end
% 
%            if p.correctcamera %apply offset and brightfield correction
%                % scmosroi 
%                if isempty(obj.offsetmapuse)
%                     loadcamcalibrationfile(obj,p,imgp);
%                end
%                imphot=(single(imgp)-obj.offsetmapuse)*obj.gainuse;   
% %                imphot=(single(imgp)-obj.offsetmap(roi(1)+1:roi(1)+roi(3),roi(2)+1:roi(2)+roi(4)))... %*obj.adu2phot...
% %                    ./obj.gainmap(roi(1)+1:roi(1)+roi(3),roi(2)+1:roi(2)+roi(4)); %XXX correct or exchanged??? a: looks fine
%            else
%                imphot=(single(imgp)-obj.offset)*obj.adu2phot;
%            end
               datao=data;
               datao.data=imphot;  
           if obj.preview && ~isempty(imphot)
               obj.setPar('preview_image',imphot);
           end
           
           %save raw imageas
           if (mod(data.frame,p.diffrawframes)==0 || data.frame==1) && ~obj.preview && p.diffrawframes>0
               if obj.rawimagecounter>length(obj.rawimagestruct)
                   obj.rawimagestruct(end+100).image=getfirstimage(imphot)*0;
                   obj.rawimagestruct(end+100).frame=-1;
               end
               obj.rawimagestruct(obj.rawimagecounter+1).image=getfirstimage(imphot);
               obj.rawimagestruct(obj.rawimagecounter+1).frame=data.frame;
               obj.rawimagecounter=obj.rawimagecounter+1;
           end
           if isempty(obj.rawaverage)
               obj.rawaverage=double(getfirstimage(imphot));
               obj.rawaveragetype=getfirstimage(imphot);
           else
               obj.rawaverage=obj.rawaverage+double(getfirstimage(imphot));
           end
           obj.rawaveragecounter=obj.rawaveragecounter+1;
        end
       
    end
end

function img= getfirstimage(imin)
s=size(imin);
if s>3
    img=imin(:,:,1,1);
else
    img=imin;
end
end

function imphot=convertsingleimage(obj,imin,p,offsetmap)
imgp=makepositive(imin);
if p.emmirror && obj.loc_cameraSettings.EMon  
        imgp=imgp(:,end:-1:1);
end
if nargin>3 && ~isempty(offsetmap)
% if p.correctcamera %apply offset and brightfield correction
   % scmosroi 
%    if isempty(obj.offsetmapuse)
%         loadcamcalibrationfile(obj,p,imgp);
%    end
   imphot=(single(imgp)-offsetmap)*obj.gainuse;   
%                imphot=(single(imgp)-obj.offsetmap(roi(1)+1:roi(1)+roi(3),roi(2)+1:roi(2)+roi(4)))... %*obj.adu2phot...
%                    ./obj.gainmap(roi(1)+1:roi(1)+roi(3),roi(2)+1:roi(2)+roi(4)); %XXX correct or exchanged??? a: looks fine
else
   imphot=(single(imgp)-obj.offset)*obj.adu2phot;
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
description=obj.loc_cameraSettings.description;
fi=obj.getPar('loc_fileinfo');
obj.loc_cameraSettings=copyfields(obj.loc_cameraSettings,fi);
for k=length(fn):-1:1
    fields{k}=description.(fn{k});
%     fields{k}=fn{k};
    if myisfield(fi,fn{k})
        defAns{k}=converttostring(fi.(fn{k}));
    else
        defAns{k}=converttostring(obj.loc_cameraSettings.(fn{k}));
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

function loadcamcalibrationfile(obj,p,imgp)
    obj.setPar('cam_varmap',[]);
    [gainmap,offsetmap,varmap,roi]=makegainoffsetCMOS(obj.loc_cameraSettings.correctionfile,obj.loc_cameraSettings.exposure);
    
    if ~isempty(gainmap)
%         obj.setPar('cam_varmap',single(varmap));
%         obj.setPar('scmos_roi',roi);
        obj.varmap=single(varmap);
        obj.scmosroi=roi;
        obj.gainmap=single(gainmap);
        obj.offsetmap=single(offsetmap);
        
        
       roiimg=p.loc_cameraSettings.roi;
       scmosroi=obj.scmosroi;

       if ~isempty(scmosroi) % specified: use 
           roi(3:4)=roiimg(3:4);
           roi(1:2)=roiimg(1:2)-scmosroi(1:2);
       elseif all(size(obj.offsetmap)==size(imgp)) % not specified: same size as image use like that
           roi=roiimg;
           roi(1:2)=0;
           disp('no scmos ROI specified but scmos calibration size equal to image size: assume it is the same ROI');
       else 
           roi=roiimg;
           roi(1:2)=roi(1:2); %zero based;
           disp('no scmos ROI specified: assume entire chip used for calibration');
       end
       sg=size(obj.gainmap);
       if any(sg([2 1])<roi(1:2)+roi(3:4)) %gainmap too small
           disp('scmos calibration map not compatible with image size.')
           roi(1:2)=0;
       end
%        roi(1)=roi(1)+5 %test
%        gainhere=(obj.gainmap(roi(1)+1:roi(1)+roi(3),roi(2)+1:roi(2)+roi(4)));
%        obj.offsetmapuse=obj.offsetmap(roi(1)+1:roi(1)+roi(3),roi(2)+1:roi(2)+roi(4));
%        obj.gainuse=median(gainhere(:));
%        varmap=obj.varmap(roi(1)+1:roi(1)+roi(3),roi(2)+1:roi(2)+roi(4));
        
%        This seems to work together with the calibrateCMOS.m plugin
       gainhere=(obj.gainmap(roi(2)+1:roi(2)+roi(4),roi(1)+1:roi(1)+roi(3)));
       obj.offsetmapuse=obj.offsetmap(roi(2)+1:roi(2)+roi(4),roi(1)+1:roi(1)+roi(3));
       obj.gainuse=median(gainhere(:));
       varmap=obj.varmap(roi(2)+1:roi(2)+roi(4),roi(1)+1:roi(1)+roi(3));
       
        obj.setPar('cam_varmap',varmap);
    
    else
        if obj.getSingleGuiParameter('correctcamera')
            disp(['could not find camera correction file ' camfile])
        end
        obj.setPar('cam_varmap',[]);
    end
   
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

function mirrorem_callback(a,b,obj)
% if ~isempty(obj.imloader)
%      fileinf=obj.imloader.metadata;
            if obj.loc_cameraSettings.EMon && obj.getSingleGuiParameter('emmirror')  %if em gain on and mirrorem on: switch roi
                obj.loc_cameraSettings.EMmirror=true;
            else 
                obj.loc_cameraSettings.EMmirror=false;
            end
%             obj.mirrorem=fileinf.EMmirror;
%             if obj.getSingleGuiParameter('padedges') 
%                 dr=obj.getSingleGuiParameter('padedgesdr');
%                 fileinf.roi(1:2)=fileinf.roi(1:2)-dr;
%                 fileinf.roi(3:4)=fileinf.roi(3:4)+2*dr;
%                 obj.edgesize=dr;
%             end
%             obj.setPar('loc_fileinfo',fileinf);
% end
end

function out=converttostring(in)
if iscell(in)
    out=join(in,',');
    out=out{1};
elseif ischar(in)
    out=in;
else
    out=num2str(in);
end

end

function pard=guidef(obj)
pard.text.object=struct('Style','text','String','Metadata:');
pard.text.position=[1,1];
pard.text.Width=.7;
pard.text.Optional=true;
 
pard.loadmetadata.object=struct('Style','pushbutton','String','Load');
pard.loadmetadata.position=[1,4.];
pard.loadmetadata.TooltipString=sprintf('Load camera settings metadata from image or _sml.mat file.');
pard.loadmetadata.Optional=true;
pard.loadmetadata.Width=0.4;
pard.calibrate.object=struct('Style','pushbutton','String','calibrate');
pard.calibrate.position=[1,4.4];
pard.calibrate.TooltipString=sprintf('calibrate gain and offset from images');
pard.calibrate.Optional=true;
pard.calibrate.Width=0.6;

pard.camparbutton.object=struct('Style','pushbutton','String','set Cam Parameters');
pard.camparbutton.position=[1,1.7];
pard.camparbutton.Width=1.3;
pard.camparbutton.TooltipString=sprintf('Edit camera acquisition parameters.');
pard.lockcampar.object=struct('Style','checkbox','String','Lock','Value',0);
pard.lockcampar.position=[1,3];
pard.lockcampar.Width=0.5;
pard.lockcampar.TooltipString=sprintf('Do not overwrite camera parameters automatically, but keep those set manually.');
pard.lockcampar.Optional=true;

pard.correctcamera.object=struct('Style','checkbox','String','Correct offset','Value',0);
pard.correctcamera.position=[2,1];
pard.correctcamera.Width=2;
pard.correctcamera.TooltipString=sprintf('Apply darkfield and brightfield correction.');
pard.correctcamera.Optional=true;

pard.emmirror.object=struct('Style','checkbox','String','Mirror if EM mode','Value',1,'Callback',{{@mirrorem_callback,obj}});
pard.emmirror.position=[2,3];
pard.emmirror.Width=2;
pard.emmirror.TooltipString=sprintf('Mirror the image if acquired with EM gain.');
pard.emmirror.Optional=true;

pard.plugininfo.type='WorkflowModule'; 
pard.plugininfo.description='Interprets metadata and converts camera ADUs into photons- Metadata can be overwritten manually or loaded from a SMAP _sml.mat data file.';
end