classdef iterativeMDfitter<interfaces.WorkflowModule
%     Single-emiiter fitter that considers all ther candidate
%     localizations. More robust for high-density localization than a
%     single-emitter fitter
    properties
        splinePSF
        splinenorm
        preview
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
            psf.roisize=p.roisize;
            obj.setPar('loc_ROIsize',p.roisize)
            psf.loadmodel(p.cal_3Dfile);
                 
            obj.splinePSF=psf;
            pim=psf.PSF(struct('x',0,'y',0,'z',0));
            obj.splinenorm=sum(pim(:))/max(pim(:));
            obj.preview=obj.getPar('loc_preview');
        end
        function outputdat=run(obj,data,p)

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
%             maxima.N=image(sub2ind(size(image),maxima.y,maxima.x))*obj.splinenorm;
%             maxima.z=0*maxima.N;
%             maxima.z=maxima.znm;
%             maxima.N=maxima.intensity;
            dr=round((obj.splinePSF.roisize-1)/2);
            maxima.z=-maxima.znm; %XXXXX 

            pos.x=max(dr+1,min(size(image,2)-dr,round(maxima.x)));
            pos.y=max(dr+1,min(size(image,1)-dr,round(maxima.y)));
            maximainit=maxima;
            maximainit.x=pos.x;
            maximainit.y=pos.y;
            
            v0=0*maxima.N;
            coord=[maxima.x-pos.x,maxima.y-pos.y,maxima.z,maxima.N,v0];
            coord0=coord;
            M=single(obj.splinePSF.render(maxima,[0 size(image,2)],[0 size(image,1)]));
%             Mstart=M;
            Mi=single(obj.splinePSF.PSF(coord));      
            [~,inds]=sort(maxima.N,'descend');
            LL=ones(length(inds),1,'single');iterations=ones(length(inds),1,'single');crlb=ones(length(inds),1,5,'single');chi2=ones(length(inds),1,'single');
            for iter=1:p.iterations
                for k=1:length(inds)
                    ih=inds(k);
                    roiim=image(pos.y(ih)-dr:pos.y(ih)+dr,pos.x(ih)-dr:pos.x(ih)+dr);
                    %update image
                    M(pos.y(ih)-dr:pos.y(ih)+dr,pos.x(ih)-dr:pos.x(ih)+dr)=M(pos.y(ih)-dr:pos.y(ih)+dr,pos.x(ih)-dr:pos.x(ih)+dr)-Mi(:,:,ih);
                    roiM=M(pos.y(ih)-dr:pos.y(ih)+dr,pos.x(ih)-dr:pos.x(ih)+dr);
                    startpar=coord(ih,:);
                    [coordf,crlb(ih,:), LL(ih), iterations(ih)]=fitsingleMD(roiim,roiM,startpar,obj.splinePSF.modelpar,p.iterationsf,p.useSEfitter);
                    coord(ih,:)=coordf;
%                     bg=coordf(5);
                    coordf(5)=0; %take out fitted bg
                    Mi(:,:,ih)=obj.splinePSF.PSF(coordf);
                    M(pos.y(ih)-dr:pos.y(ih)+dr,pos.x(ih)-dr:pos.x(ih)+dr)=M(pos.y(ih)-dr:pos.y(ih)+dr,pos.x(ih)-dr:pos.x(ih)+dr)+Mi(:,:,ih);
                    dim=(roiM+Mi(:,:,ih)-roiim+coord(ih,5)).^2./(roiM+Mi(:,:,ih)+coord(ih,5));
                    chi2(k)=(mean(dim(:)));
                end   
            end
            locout=coord2loc(coord,crlb,LL,iterations,pos,data{1}.frame,chi2);
%             locout=copyfields(maxima,locout);
            locout=copyfields(locout,maxima,{'xcnn','ycnn','zcnn','prob','dx','dy'});
            if obj.preview
                
                 figure(87);
                 subplot(2,2,2)
                 imagesc(M-image); colorbar
                 subplot(2,2,1)
%                  figure(90)
                 hold off
                 imagesc(image);colorbar
                 hold on
                 plot(coord0(:,1)+pos.x,coord0(:,2)+pos.y,'mo');
                 plot(coord(:,1)+pos.x,coord(:,2)+pos.y,'k+');%,maxima.x,maxima.y,'md')
                gt=obj.getPar('loc_gt_preview');
                if ~isempty(gt)
                    plot(gt.x,gt.y,'k+')
%                     subplot(2,2,3)
%                     plot(gt.x,gt.z,'k+',coord(:,1)+pos.x,coord(:,3),'ro')
%                     xlim([0,size(image,2)])
                end
                axis('equal')
                colormap parula
                figure(87)
                subplot(2,2,4)
                plot(coord(:,1)+pos.x,-LL/((2*dr+1)^2),'k+',coord(:,1)+pos.x,chi2,'x')
                xlim([0,size(image,2)])
    %             pause(.5)

%                 locsfit.x=coord(:,1)+pos.x; locsfit.y=coord(:,2)+pos.y;locsfit.z=coord(:,3);locsfit.N=coord(:,4);locsfit.bg=coord(:,5)*0;
%                 gt.bg=0; 
            end
            outputdat=data{1};
            outputdat.data=locout;
        end
    end
end

function locout=coord2loc(coord,crlb,LL,iterations,pos,frame,chi2)
locout.xpix=single(coord(:,1)+pos.x);
locout.ypix=single(coord(:,2)+pos.y);
locout.znm=-single(coord(:,3)); %XXXXX
locout.phot=single(coord(:,4));
locout.bg=single(coord(:,5));
locout.xpixerr=single(sqrt(crlb(:,1)));
locout.ypixerr=single(sqrt(crlb(:,2)));
locout.zerr=single(sqrt(crlb(:,3)));
locout.photerr=single(sqrt(crlb(:,4)));
locout.bgerr=single(sqrt(crlb(:,5)));
locout.logLikelihood=single(LL);
locout.iterations=single(iterations);
locout.frame=zeros(size(locout.xpix))+frame;
locout.chi2=single(chi2);

end
function [coord,crlb, LogL, iterations]=fitsingleMD(roiim,roiM,startpar,spline,iterf,useSEfitter)
dr=round((size(roiim,1)-1)/2);

%    cor(:,3)=-locs(:,3)/obj.modelpar.dz+obj.modelpar.z0;
initp(3:4)=startpar(4:5);
initp(5)=-startpar(3)/spline.dz+spline.z0;
initp(1:2)=startpar([2 1])+dr;
if useSEfitter
    [P,CRLB,LogL]=mleFit_LM(roiim,5,50,spline.coeff);
else
 [P,CRLB,LogL]=mleFit_LM_HD_SE(roiim,iterf,spline.coeff,roiM,initp);
end
%   
 coord(4:5)=P(3:4);
coord(1:2)=P([2 1])-dr;
coord(3)=-(P(5)-spline.z0)*spline.dz;

crlb=CRLB([2 1 5 3 4]);
crlb(3)=crlb(3)*spline.dz^2;
iterations=P(end);

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

pard.roisizet.object=struct('Style','text','String','Roi size');
pard.roisizet.position=[3,1];
pard.roisizet.Width=1;

pard.roisize.object=struct('Style','edit','String','13');
pard.roisize.position=[3,2];
pard.roisize.Width=0.5;

pard.iterationst.object=struct('Style','text','String','Iterations repeat');
pard.iterationst.position=[3,3];
pard.iterationst.Width=1;

pard.iterations.object=struct('Style','edit','String','3');
pard.iterations.position=[3,4];
pard.iterations.Width=0.5;

pard.iterationsft.object=struct('Style','text','String','Iterations fit');
pard.iterationsft.position=[4,1];
pard.iterationsft.Width=1;

pard.iterationsf.object=struct('Style','edit','String','30');
pard.iterationsf.position=[4,2];
pard.iterationsf.Width=0.5;

pard.useSEfitter.object=struct('Style','checkbox','String','use SE fitter','Value',0);
pard.useSEfitter.position=[5,1];
pard.useSEfitter.Width=1;

pard.syncParameters={{'cal_3Dfile','cal_3Dfile',{'String'}}};

pard.plugininfo.type='WorkflowModule'; 
pard.plugininfo.description='Single-emiiter fitter that considers all ther candidate localizations. More robust for high-density localization than a single-emitter fitter';
end