classdef iterativeMDfitter<interfaces.WorkflowModule
    properties
        splinePSF
        splinenorm
    end
    methods
        function obj=iterativeMDfitter(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=2; 
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
            obj.setInputChannels(2,'frame');
        end
        function prerun(obj,p)  
            p=obj.getAllParameters;
      
            psf=splinePSF;
            psf.loadmodel(p.cal_3Dfile);
                 
            obj.splinePSF=psf;
            pim=psf.PSF(struct('x',0,'y',0,'z',0));
            obj.splinenorm=sum(pim(:))/max(pim(:));
        end
        function outputdat=run(obj,data,p)
            outputdat=[];
            image=data{1}.data;%get;
            if isempty(image) 
                outputdat=data{1};
                return
            end
            maxima=data{2}.data;%get;
            if isempty(maxima.x)
                dato=data{1};%.copy;
                dato.data=[];
                outputdat=dato;
                return;
            end
            maxima.N=image(sub2ind(size(image),maxima.y,maxima.x))*obj.splinenorm;
            maxima.z=0*maxima.N;
            dr=round((obj.splinePSF.roisize-1)/2);

            pos.x=max(1,min(size(image,2)-dr,round(maxima.x)));
            pos.y=max(1,min(size(image,1)-dr,round(maxima.y)));
            
            v0=0*maxima.N;
            coord=[pos.x-maxima.x,pos.y-maxima.y,v0,maxima.N,v0];
            M=obj.splinePSF.render(maxima,[0 size(image,1)],[0 size(image,2)]);
            Mi=obj.splinePSF.PSF(coord);      
            [~,inds]=sort(maxima.N,'descend');
            for k=1:length(inds)
                ih=inds(k);
                roiim=image(pos.y(ih)-dr:pos.y(ih)+dr,pos.x(ih)-dr:pos.x(ih)+dr);
                %update image
                M(pos.y(ih)-dr:pos.y(ih)+dr,pos.x(ih)-dr:pos.x(ih)+dr)=M(pos.y(ih)-dr:pos.y(ih)+dr,pos.x(ih)-dr:pos.x(ih)+dr)-Mi(:,:,ih);
                roiM=M(pos.y(ih)-dr:pos.y(ih)+dr,pos.x(ih)-dr:pos.x(ih)+dr);
                startpar=coord(ih,:);
                coordf=fitsingleMD(roiim,roiM,startpar);
                Mi(:,:,ih)=obj.splinePSF.render(coordf,[0 size(image,1)],[0 size(image,2)]);
                M(pos.y(ih)-dr:pos.y(ih)+dr,pos.x(ih)-dr:pos.x(ih)+dr)=M(pos.y(ih)-dr:pos.y(ih)+dr,pos.x(ih)-dr:pos.x(ih)+dr)+Mi(:,:,ih);
            end       
        end
        

    end
end

function loadcall_callback(a,b,obj)
p=obj.getAllParameters;
if isempty(p.cal_3Dfile)
    path=obj.getGlobalSetting('DataDirectory');
    fh=obj.getPar('loc_fileinfo');
    if ~isempty(fh)
        path=fileparts(fh.imagefile);
    end  
    p.cal_3Dfile=[path filesep '*3dcal.mat'];
end
[f,p]=uigetfile(p.cal_3Dfile);
if f
    l=load([p f]);
    if ~isfield(l,'outforfit') && ~isfield(l,'SXY') && ~isfield(l,'cspline')
        msgbox('no 3D data recognized. Select other file.');
    end
    obj.setGuiParameters(struct('cal_3Dfile',[p f]));
    
end
end

function pard=guidef(obj)
pard.loadcal.object=struct('Style','pushbutton','String','Load 3D cal','Callback',{{@loadcall_callback,obj}});
pard.loadcal.position=[2,1];
pard.loadcal.Width=.75;
pard.cal_3Dfile.object=struct('Style','edit','String','');
pard.cal_3Dfile.position=[2,1.75];
pard.cal_3Dfile.Width=1.75;
pard.cal_3Dfile.TooltipString=sprintf('3D calibration file for astigmtic 3D. \n Generate from bead stacks with plugin: Analyze/sr3D/CalibrateAstig');


% pard.syncParameters={{'loc_ROIsize','loc_ROIsize',{'String'}},{'loc_filterforfit','loc_filterforfit',{'String'}}};

pard.plugininfo.type='WorkflowModule'; 
pard.plugininfo.description='This plugin cuts out regions of interest of a defined size around the candidate positions and passes these on to the fitter';
end