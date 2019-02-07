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
            psf.roisize=17;
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

            pos.x=max(dr+1,min(size(image,1)-dr,round(maxima.x)));
            pos.y=max(dr+1,min(size(image,2)-dr,round(maxima.y)));
            maximainit=maxima;
            maximainit.x=pos.x;
            maximainit.y=pos.y;
            
            v0=0*maxima.N;
            coord=[maxima.x-pos.x,maxima.y-pos.y,v0,maxima.N,v0];
            M=single(obj.splinePSF.render(maxima,[0 size(image,2)],[0 size(image,1)]));
            Mstart=M;
            Mi=single(obj.splinePSF.PSF(coord));      
            [~,inds]=sort(maxima.N,'descend');
            for iter=1:5
            for k=1:length(inds)
                ih=inds(k);
                roiim=image(pos.y(ih)-dr:pos.y(ih)+dr,pos.x(ih)-dr:pos.x(ih)+dr);
                %update image
                M(pos.y(ih)-dr:pos.y(ih)+dr,pos.x(ih)-dr:pos.x(ih)+dr)=M(pos.y(ih)-dr:pos.y(ih)+dr,pos.x(ih)-dr:pos.x(ih)+dr)-Mi(:,:,ih);
                roiM=M(pos.y(ih)-dr:pos.y(ih)+dr,pos.x(ih)-dr:pos.x(ih)+dr);
                startpar=coord(ih,:);
                coordf=fitsingleMD(roiim,roiM,startpar,obj.splinePSF.modelpar);
                coord(ih,:)=coordf;
                bg=coordf(5);
                coordf(5)=0; %take out fitted bg
                Mi(:,:,ih)=obj.splinePSF.PSF(coordf);
                M(pos.y(ih)-dr:pos.y(ih)+dr,pos.x(ih)-dr:pos.x(ih)+dr)=M(pos.y(ih)-dr:pos.y(ih)+dr,pos.x(ih)-dr:pos.x(ih)+dr)+Mi(:,:,ih);
            end   
            end
             figure(87);
             subplot(2,2,2)
             imagesc(M-image); colorbar
             subplot(2,2,1)
             hold off
             imagesc(image);
             hold on
            gt=obj.getPar('loc_gt_preview');
            plot(gt.x,gt.y,'ko',coord(:,1)+pos.x,coord(:,2)+pos.y,'kx',maxima.x,maxima.y,'md')
            subplot(2,2,3)
            plot(gt.x,gt.z,'ko',coord(:,1)+pos.x,coord(:,3),'kx')
            xlim([0,size(image,2)])
%             pause(.5)
            
            locsfit.x=coord(:,1)+pos.x; locsfit.y=coord(:,2)+pos.y;locsfit.z=coord(:,3);locsfit.N=coord(:,4);locsfit.bg=coord(:,5)*0;
            gt.bg=0; 
%             Mgt=single(obj.splinePSF.render(gt,[0 size(image,2)],[0 size(image,1)]));
%              subplot(2,2,4);
%              imagesc(Mgt-image); colorbar
%              drawnow
%              pause(1)
            
          
        end
        

    end
end

function coord=fitsingleMD(roiim,roiM,startpar,spline)
dr=round((size(roiim,1)-1)/2);

%    cor(:,3)=-locs(:,3)/obj.modelpar.dz+obj.modelpar.z0;
initp(3:4)=startpar(4:5);
initp(5)=-startpar(3)/spline.dz+spline.z0;
initp(1:2)=startpar([2 1])+dr;
 [P,CRLB,LogL]=mleFit_LM_HD_SE(roiim,50,spline.coeff,roiM,initp);
%   [P,CRLB,LogL]=mleFit_LM(roiim,5,50,spline.coeff);
 coord(4:5)=P(3:4);
coord(1:2)=P([2 1])-dr;
coord(3)=-(P(5)-spline.z0)*spline.dz;


end

function loadcall_callback(a,b,obj)
p=obj.getAllParameters;
if isempty(p.cal_3Dfile)
    path=obj.getGlobalSetting('DataDirectory');
    fh=obj.getPar('loc_fileinfo');
    if ~isempty(fh) && ~isempty(fh.imagefile)
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


pard.syncParameters={{'cal_3Dfile','cal_3Dfile',{'String'}}};

pard.plugininfo.type='WorkflowModule'; 
pard.plugininfo.description='This plugin cuts out regions of interest of a defined size around the candidate positions and passes these on to the fitter';
end