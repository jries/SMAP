classdef GetIntensitiesSALM<interfaces.DialogProcessor
    % gets intensities from camera images at positions of localizations and
    % at transformed positions
    properties (Access=private)
        figure
    end
    methods
        function obj=GetIntensitiesSALM(varargin)   
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.history=true;
            obj.showresults=false;
        end
        
        function out=run(obj,p)
            out=[];
            % make T
            pt.uselayers=true;
            pt.reflayer=p.ualayer; pt.targetlayer=p.salayer;
            pt.currentfileinfo=obj.locData.files.file(1).info;
            pt.useT=true;
            cal3D=load(p.cal_3Dfile);
            transform=cal3D.transformation;
            pt.Tfile=transform;
            pt.register_parameters.pixelsizenm=200;
            pt.register_parameters.maxshift_corr=10000;
            pt.register_parameters.maxshift_match=350;
            pt.resultstabgroup=obj.resultstabgroup;
            pt.register_parameters.maxlocsused=50000;
            if p.makeT
            transform=transform_locsN(obj.locData,pt);
            end
%             [locsa,xx]=obj.locData.getloc({'xnm','ynm','frame','filenumber','znm','PSFxnm'},'layer',p.salayer.Value,'position','roi','grouping','ungrouped');
%             [locua,yy]=obj.locData.getloc({'xnm','ynm','frame','filenumber','znm','PSFxnm'},'layer',p.ualayer.Value,'position','roi','grouping','ungrouped');
%             pos=transformcal.transformToReference(2,horzcat(locua.xnm,locua.ynm),'nm');
%     
            %make WF
            if isempty(obj.figure)||~isvalid(obj.figure)
                obj.figure=figure;
            end
            f=obj.figure;
            f.Visible='on';
%             wffile='settings/workflows/get2CIntensityImagesWF_group.mat';
%             wffile='settings/workflows/get2CIntensityImagesWF2';
            wffile='settings/workflows/get2CIntensityImagesWF3';
%             wffile='settings/workflows/get2CIntensityImagesWF_group_bg';
            wf=interfaces.Workflow(f,obj.P);
            wf.attachLocData(obj.locData);
            wf.makeGui;
            wf.load(wffile);
            
      
            p.loaderclass=wf.module('TifLoader').getSingleGuiParameter('loaderclass'); 
            p.loaderclass.Value=1;%single MM
            p.framestop=max(obj.locData.loc.frame);
            wf.module('TifLoader').setGuiParameters(p);
                  
            
            file=obj.locData.files.file;
            wf.module('IntLoc2posN').filestruc=file;
            wf.module('IntLoc2posN').transformation=transform;
            wf.module('IntLoc2posN').setGuiParameters(p);
            
            pe=obj.children.evaluate.getGuiParameters(true);
            bgmode=pe.children.panel_3.fitbg.Value;
            p.calculatebg=bgmode==3; 
            wf.module('Roi_bg').setGuiParameters(p);          

            
            wf.module('EvaluateIntensity_s').setGuiParameters(pe,true);
            obj.setPar('loc_preview',false);
            
            evalfit=wf.module('EvaluateIntensity_s').children.panel_3;
            rsfit=evalfit.getSingleGuiParameter('roisize_fit');
            pfit=evalfit.getAllParameters;
              pfit.fixz0=1;
             evalfit.setGuiParameters(pfit)
            
            
            p.loc_ROIsize=rsfit+2;
%             p.loc_fitgrouped=p.fitgrouped;
%             wf.module('RoiCutterWF_groupExt').setGuiParameters(p);
            wf.module('RoiCutterWF').setGuiParameters(p);

            
            
            % now first SA
    
                obj.setPar('intensity_channel','s')
                wf.module('TifLoader').addFile(p.tiffiletarget,true);   
                wf.module('TifLoader').setGuiParameters(struct('mirrorem',p.mirrorem))
%                 wf.module('EvaluateIntensity_s').extension='t';

                overwritedefaultcamera(obj)
                wf.module('IntLoc2posN').setGuiParameters(struct('transformtotarget',true));
                wf.run;
        %then fit UA if needed with modified settings
            if p.evalua
                obj.setPar('intensity_channel','u')
%                 if isempty(p.tiffileref) %when same, dont select again
                    p.tiffileref=p.tiffiletarget;
                    p.mirroremref=p.mirrorem;
%                 end
                
                pfit.fixz0=0;
                evalfit.setGuiParameters(pfit)
             
                wf.module('TifLoader').addFile(p.tiffileref,true);   
                wf.module('TifLoader').setGuiParameters(struct('mirrorem',p.mirroremref))
%                 wf.module('EvaluateIntensity_s').extension='r';
                
                wf.module('IntLoc2posN').setGuiParameters(struct('transformtotarget',false));
                overwritedefaultcamera(obj)
                wf.run;
            end
            delete(f);
            fo=strrep(obj.locData.files.file(1).name,'_sml.mat','_dc_sml.mat');
            obj.locData.savelocs(fo);
            obj.locData.regroup;
            obj.setPar('locFields',fieldnames(obj.locData.loc))

        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function initGui(obj)
            par.Vpos=3;
            par.Xpos=3;
            obj.children.evaluate=makeplugin(obj,{'WorkflowModules','IntensityCalculator','EvaluateIntensity_s'},obj.handle,par);
%             obj.guihandles.loadbutton.Callback=@obj.loadbutton;

        end
        function loadbutton_3dcal(obj,a,b)
            file3dtest=obj.guihandles.cal_3Dfile.String;
            if isempty(file3dtest)
                file3dtest=obj.locData.files.file(1).savefit.fitparameters.MLE_GPU_Yiming.cal_3Dfile;
                if ~exist(file3dtest,'file')
                    file3dtest=getrawtifpath(obj.locData);
                    file3dtest=[fileparts(file3dtest) filesep '*_3dcal.mat'];
                end
            end
            [f,path]=uigetfile(file3dtest,'Select 3D calibration');
            if f
                obj.guihandles.cal_3Dfile.String=[path f];
                obj.setPar('cal_3Dfile',[path f])
            end      
        end
        function loadbutton_tif(obj,a,b,field)
            fn= obj.guihandles.(field).String;
            if isempty(fn)
                fn=getrawtifpath(obj.locData);
            end
            
            [f,path]=uigetfile(fn,'Select raw tiff file');
            if f
                obj.guihandles.(field).String=[path f];
            end      
            r=imageloaderAll([path f],[],obj.P);
            mirror=r.metadata.EMon;
            r.close;
            obj.setGuiParameters(struct('mirrorem',mirror));
            
        end
    end
end




function pard=guidef(obj)


pard.texttr.object=struct('String','SA:','Style','text');
pard.texttr.position=[1,2.5];
pard.texttr.Width=0.25;

pard.salayer.object=struct('Style','popupmenu','String','layer1|layer2|layer3|layer4|');
pard.salayer.position=[1,2.65];
pard.salayer.object.TooltipString='layer';
pard.salayer.Width=0.85;

pard.texttl.object=struct('String','UA:','Style','text');
pard.texttl.position=[1,3.5];
pard.texttl.Width=0.25;

pard.ualayer.object=struct('Style','popupmenu','String','layer1|layer2|layer3|layer4|','Value',1);
pard.ualayer.position=[1,3.65];
pard.ualayer.object.TooltipString='layer';
pard.ualayer.Width=0.85;
ol={'texttr','salayer','texttl','ualayer'};
p(1).value=0; p(1).on={}; p(1).off=ol;
p(2).value=1; p(2).on=ol; p(2).off={};

pard.makeT.object=struct('Style','checkbox','String','make Transform','Value',1,'Callback',{{@obj.switchvisible,p}});
pard.makeT.position=[1,1];
pard.makeT.Width=1.5;
  

pard.tiffiletarget.object=struct('Style','edit','String','');
pard.tiffiletarget.position=[2,1.6];
pard.tiffiletarget.Width=2.7;

pard.loadbuttontiftarget.object=struct('Style','pushbutton','String','load tif','Callback',{{@obj.loadbutton_tif,'tiffiletarget'}});
pard.loadbuttontiftarget.position=[2,4.3];
pard.loadbuttontiftarget.Width=0.7;

pard.mirrorem.object=struct('Style','checkbox','String','EM','Value',1);
pard.mirrorem.position=[2,1.0];
pard.mirrorem.Optional=true;
pard.mirrorem.Width=0.6;


pard.cal_3Dfile.object=struct('Style','edit','String','');
pard.cal_3Dfile.position=[3,1.6];
pard.cal_3Dfile.Width=2.7;
pard.cal_3Dfile.TooltipString=sprintf('transformation file. Created e.g. with RegisterLocs');

pard.loadbuttoncal.object=struct('Style','pushbutton','String','load 3dcal','Callback',@obj.loadbutton_3dcal);
pard.loadbuttoncal.position=[3,4.3];
pard.loadbuttoncal.Width=0.7;
% pard.loadbuttoncal.TooltipString=pard.Tfile.TooltipString;

% 
pard.evalua.object=struct('Style','checkbox','String','evaluate also UA','Value',1);
pard.evalua.position=[4,1];
pard.evalua.Width=2;

pard.loc_fitgrouped.object=struct('Style','checkbox','String','evaluate grouped images','Value',0);
pard.loc_fitgrouped.position=[4,3];
pard.loc_fitgrouped.Width=2;



p(1).value=0; p(1).on={}; p(1).off={'bgfunction','bgfunctionpar','t1','numframes_bg'};
p(2).value=1; p(2).on=p(1).off; p(2).off={};
pard.calculatebg.object=struct('Style','checkbox','String','calculate BG','Callback',{{@obj.switchvisible,p}});
pard.calculatebg.position=[5,1];
pard.calculatebg.Width=1.2;

pard.bgfunction.object=struct('Style','popupmenu','String',{{'quantile'}},'Visible','off');
pard.bgfunction.position=[5,2.2];
pard.bgfunction.Width=0.9;
pard.bgfunctionpar.object=struct('Style','edit','String','0.5','Visible','off');
pard.bgfunctionpar.position=[5,3.1];
pard.bgfunctionpar.Width=0.3;

pard.t1.object=struct('Style','text','String','# frames','Visible','off');
pard.t1.position=[5,3.8];
pard.t1.Width=.9;
pard.numframes_bg.object=struct('Style','edit','String','20','Visible','off');
pard.numframes_bg.position=[5,4.7];
pard.numframes_bg.Width=0.3;


pard.syncParameters={{'cal_3Dfile','cal_3Dfile',{'String'}}};
pard.plugininfo.type='ProcessorPlugin';
% pard.plugininfo.description=sprintf(['This plugin gets intensities from camera images at positions of localizations and at transformed positions \n',...
%     'This plugin uses a transformation to find for every localization the position in the other channel and then determines the intensity in both channels.\n',...
%     '1.	Load a transformation\n',...
%     '2.	Per default, this plugin does median filtering. Select the spatial and temporal spacing for this (dx, dt).\n',...
%     '3.	Select one or several plugins which determine the intensity:\n',...
%     '\t a.	Roi2int_sum: uses a ROI (set size) to determine intensity, and a larger ROI for the background.\n',...
%     '\t b.	Roi2int_fit: Uses a Gaussian fit to determine intensity and background. The position is fixed to the fitted position. You can use the fitted PSF size or fix it. If fit on BG is checked, the background is subtracted prior to fitting and the fit is performed with background set to zero. Otherwise the background is a fitting parameter.\n',...
%     '4.	Press Run and when asked select the original raw camera images. The results are automatically saved with the _dc in the file name.\n']);
pard.plugininfo.name='Intensities for SALM';
% pard.plugininfo.description
end



function tiffile=getrawtifpath(locData)
    if isfield(locData.files.file(1).info,'imagefile')
        tiffile=locData.files.file(1).info.imagefile;
    else
        tiffile=locData.files.file(1).info.basefile;
    end
    if ~exist(tiffile,'file')
        disp('Get2CIntImagesWF ine 40: check if it works')
        tiffile=strrep(tiffile,'\','/');
        ind=strfind(tiffile,'/');
        for k=1:length(ind)
            tiffileh=[tiffile(1:ind(k)) '_b_' tiffile(ind(k)+1:end)];
            if exist(tiffileh,'file')
                tiffile=tiffileh;
            end
        end
    end
    if ~exist(tiffile,'file')
        tiffile=strrep(locData.files.file(1).name,'_sml.mat','.tif');
    end
end

function plugino=makeplugin(obj,pluginpath,handle,pgui)
plugino=plugin(pluginpath{:});
plugino.attachPar(obj.P);
plugino.setGuiAppearence(pgui);
plugino.attachLocData(obj.locData);
plugino.handle=handle;
%XXXX hack
plugino.makeevaluators(1:2)=false;
plugino.makeGui;
pluginname=pluginpath{end};
obj.children.(pluginname)=plugino;
end

function overwritedefaultcamera(obj)
loc_cameraSettings=obj.getPar('loc_cameraSettings');
if strcmp(loc_cameraSettings.comment,'Default camera')
      disp('Default camera recognized. Overwrite settings with settings from locData')
    loc_cameraSettingsnew=copyfields([],obj.locData.files.file(1).info,fieldnames(loc_cameraSettings));
    obj.setPar('loc_fileinfo',loc_cameraSettingsnew);
    obj.setPar('loc_fileinfo_set',loc_cameraSettingsnew);
end
end