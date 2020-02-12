classdef iterativeMDfitter<interfaces.WorkflowModule
%     Single-emiiter fitter that considers all ther candidate
%     localizations. More robust for high-density localization than a
%     single-emitter fitter
    properties
        splinePSF
        splinenorm
        preview
        EMon
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
%             psf.roisize=p.roisize;
            obj.setPar('loc_ROIsize',p.roisize)
            psf.loadmodel(p.cal_3Dfile);
            psf.roisize=min(max(15,round(1.3*p.roisize*2)+1),  size(psf.modelpar.coeff,1)-3);
                 
            obj.splinePSF=psf;
            pim=psf.PSF(struct('x',0,'y',0,'z',0));
            obj.splinenorm=sum(pim(:))/max(pim(:));
            obj.preview=obj.getPar('loc_preview');
            cs=obj.getPar('loc_cameraSettings');
            obj.EMon=cs.EMon+1;
        end
        function outputdat=run(obj,data,p)

            image=data{1}.data;%get;
            if isempty(image) 
                outputdat=data{1};
                return
            end
            maxima=data{2}.data;%get;
            if isempty(maxima.xpix)
                dato=data{1};%.copy;
                dato.data=[];
                outputdat=dato;
                return;
            end
            sim=size(image);
            
            %render image and errors with maximum PSF size!
            
%             maxima.N=image(sub2ind(size(image),maxima.ypix,maxima.xpix))*obj.splinenorm;
%             maxima.z=0*maxima.N;
%             maxima.z=maxima.znm;
%             maxima.N=maxima.intensity;
            drfit=round((p.roisize-1)/2);
            drpsf=round((obj.splinePSF.roisize-1)/2);
            spsf=obj.splinePSF.roisize;
            ddr=drpsf-drfit;
%             calculatebg=true;
%             bgestimate=min(image(:));
            sizepixfit=(2*drfit+1)^2;
            

            sigma2=1.5^2; %estimated sigma of PSF squared
            pos.x=max(drfit+1,min(size(image,2)-drfit,round(maxima.xpix)));
            pos.y=max(drfit+1,min(size(image,1)-drfit,round(maxima.ypix)));
            maximainit=maxima; %what are these needed for???
            maximainit.x=pos.x;
            maximainit.y=pos.y;  
            if p.estimateStartParameters
                v0=0*maxima.xpix;
                bgglobal=quantile(image(:),0.5);
                maxima.znm=v0;
                maxima.bg=v0;
                maxima.phot=v0;
                for k=1:length(maxima.xpix)
                    roiim=image(pos.y(k)-drfit:pos.y(k)+drfit,pos.x(k)-drfit:pos.x(k)+drfit);
%                     maxima.bg(k)=quantile(roiim(:),0.2);
                    maxima.bg(k)=bgglobal;
%                     maxima.N(k)=(sum(roiim(:))-maxima.bg(k)*sizepixfit);
%                     maxima.phot(k)=(max(roiim(:))-maxima.bg(k)*sigma2*pi);
                     maxima.phot(k)=(max(roiim(:))-maxima.bg(k))*sigma2*pi;
                end

            else %use this now for deepStorm
                maxima.znm=-maxima.znm; %XXXXX 
%                 maxima.bg=bgglobal; %currently not calculated
            end
            
            if obj.preview
                previewcollage=zeros(3*(2*drfit+1),2*(2*drfit+1),length(maxima.xpix),'single');
            end
           
            coord=[maxima.xpix-pos.x,maxima.ypix-pos.y,maxima.znm,maxima.phot,maxima.bg]; %BG set to zero
            % in iterations, previous bg used for starting value. BG always
            % as offset of single molecule model.
            coord0=coord;
            maximabg0=maxima;maximabg0.bg=maximabg0.bg*0;
            M=single(obj.splinePSF.render(maximabg0,[0 size(image,2)],[0 size(image,1)])); %start image without background. 
%             Mstart=M;
            coordbg0=coord;coordbg0(:,5)=0;
            Mi=single(obj.splinePSF.PSF(coordbg0));      
            [~,inds]=sort(maxima.phot,'descend');
            LL=ones(length(inds),1,'single');iterations=ones(length(inds),1,'single');crlb=ones(length(inds),1,5,'single');chi2=ones(length(inds),1,'single');
            for iter=1:p.iterations
                for k=1:length(inds)
                    if iter==1 && p.estimateStartParameters
                        drfith=round((p.roistart-1)/2);
                        ddrh=drpsf-drfith;
                        dzstart=single(p.zstart);
                    else
                        dzstart=0;
                        drfith=drfit;
                        ddrh=ddr;
                    end
                    ih=inds(k);
                    roiim=image(pos.y(ih)-drfith:pos.y(ih)+drfith,pos.x(ih)-drfith:pos.x(ih)+drfith);
                    %update image
                    %calculate ranges 
                    
                    yrange=max(1,pos.y(ih)-drpsf):min(pos.y(ih)+drpsf,sim(1)); 
                    yrangepsf=yrange(1)-(pos.y(ih)-drpsf)+1:yrange(end)-(pos.y(ih)+drpsf)+spsf;
                    xrange=max(1,pos.x(ih)-drpsf):min(pos.x(ih)+drpsf,sim(2)); 
                    xrangepsf=xrange(1)-pos.x(ih)+drpsf+1:xrange(end)-pos.x(ih)-drpsf+spsf;
                    
                    M(yrange,xrange)=M(yrange,xrange)-Mi(yrangepsf,xrangepsf,ih);
                    roiM=M(pos.y(ih)-drfith:pos.y(ih)+drfith,pos.x(ih)-drfith:pos.x(ih)+drfith);
                    startpar=coord(ih,:);
%                     if iter==1 && calculatebg %if no bg is passed on, use the minimum as an estimate.
%                         startpar(5)=min(roiim(:));
%                     end
                    
                    EMexcess=obj.EMon;
%                     EMexcess=1;
                    [coordf,crlb(ih,:), LL(ih), iterations(ih)]=fitsingleMD(roiim,roiM,startpar,obj.splinePSF.modelpar,p.iterationsf,p.useSEfitter,EMexcess,dzstart);
                    coord(ih,:)=coordf;
                    bgh=coordf(5);
                    coordf(5)=0; %take out fitted bg
                    Mi(:,:,ih)=obj.splinePSF.PSF(coordf);
                    M(yrange,xrange)=M(yrange,xrange)+Mi(yrangepsf,xrangepsf,ih);
                    dim=(roiM+Mi(1+ddrh:end-ddrh,1+ddrh:end-ddrh,ih)-roiim+coord(ih,5)).^2./(roiM+Mi(1+ddrh:end-ddrh,1+ddrh:end-ddrh,ih)+coord(ih,5));
                    chi2(k)=(mean(dim(:)));
%                     disp([num2str(iter) '. ' num2str(coordf(1:2)-startpar(1:2),2) ' ,dN/N=' num2str((coordf(4)-startpar(4))/coordf(4)) ' ,dz=' num2str((coordf(3)-startpar(3)))])
                    if obj.preview
                        previewcollage(1:3*(2*drfith+1),1:2*(2*drfith+1),ih)=(vertcat(horzcat(roiim-Mi(1+ddrh:end-ddrh,1+ddrh:end-ddrh,ih),roiM+bgh),horzcat(roiim,roiM+Mi(1+ddrh:end-ddrh,1+ddrh:end-ddrh,ih)+bgh),horzcat(roiim-roiM,Mi(1+ddrh:end-ddrh,1+ddrh:end-ddrh,ih)+bgh)));
              
                    end
                    
                end   
            end
            locout=coord2loc(coord,crlb,LL,iterations,pos,data{1}.frame,chi2);
%             locout=copyfields(maxima,locout);
            locout=copyfields(locout,maxima,{'xcnn','ycnn','zcnn','prob','dx','dy','photcnn'});
            if obj.preview
                f=figure(89);
                imx(previewcollage,'Parent',f)
%                 figure(89);imagesc(vertcat(horzcat(roiim-Mi(1+ddrh:end-ddrh,1+ddrh:end-ddrh,ih),roiM+bgh),horzcat(roiim,roiM+Mi(1+ddrh:end-ddrh,1+ddrh:end-ddrh,ih)+bgh),horzcat(roiim-roiM,Mi(1+ddrh:end-ddrh,1+ddrh:end-ddrh,ih)+bgh)))
                 figure(87);
                 subplot(2,2,2)
                 medbg=median(coord(:,5));
                 imagesc(M-image+medbg); colorbar
                 title('M-image+median(bg)')
                 subplot(2,2,3)
                 imagesc(M); colorbar
                 title('M')
                 subplot(2,2,1)
%                  figure(90)
                 hold off
                 imagesc(image);colorbar
                 hold on
                 
                 plot(coord0(:,1)+pos.x,coord0(:,2)+pos.y,'mo', 'DisplayName','init');
                 plot(coord(:,1)+pos.x,coord(:,2)+pos.y,'k+', 'DisplayName','fit');%,maxima.xpix,maxima.ypix,'md')
                 legend()
                 title('image')
%                  plot(pos.x,pos.y,'ks');
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
                plot(coord(:,1)+pos.x,-LL/((2*drfit+1)^2),'k+',coord(:,1)+pos.x,chi2,'x')
                
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
function [coord,crlb, LogL, iterations]=fitsingleMD(roiim,roiM,startpar,spline,iterf,useSEfitter,EMexcess,dzstart)
dr=round((size(roiim,1)-1)/2);
%    cor(:,3)=-locs(:,3)/obj.modelpar.dz+obj.modelpar.z0;
initp(3:4)=startpar(4:5)/EMexcess;
initp(5)=-startpar(3)/spline.dz+spline.z0;
initp(1:2)=startpar([2 1])+dr;
LogLold=-inf;

if useSEfitter
    [P,CRLB,LogL]=mleFit_LM(roiim/EMexcess,5,50,spline.coeff,[],[],dzstart);
else
    for k=1:length(dzstart)
        initph=single(initp);
        initph(3)=initph(3)+dzstart(k);

        [P,CRLB,LogL]=mleFit_LM_HD_SE(roiim/EMexcess,iterf,spline.coeff,roiM/EMexcess,initph);

        if LogL>LogLold
            Pf=P;
            CRLBf=CRLB;
            LogLf=LogL;
        end
    end
    P=Pf;
    CRLB=CRLBf;
    LogL=LogLf;
end

%   
 coord(4:5)=P(3:4)*EMexcess;
coord(1:2)=P([2 1])-dr;
coord(3)=-(P(5)-spline.z0)*spline.dz;

crlb=CRLB([2 1 5 3 4]);
crlb(3)=crlb(3)*spline.dz^2;
crlb(4:5)=crlb(4:5)*EMexcess;
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
    obj.setPar('cal_3Dfile',[p f]);
    
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

pard.estimateStartParameters.object=struct('Style','checkbox','String','estimate start par','Value',0);
pard.estimateStartParameters.position=[6,1];
pard.estimateStartParameters.Width=1.5;

pard.zstartt.object=struct('Style','text','String','1st zstart');
pard.zstartt.position=[6,2.5];
pard.zstartt.Width=.5;
pard.zstart.object=struct('Style','edit','String','-300 0 300');
pard.zstart.position=[6,3];
pard.zstart.Width=1;

pard.roistartt.object=struct('Style','text','String','1st ROI');
pard.roistartt.position=[6,4];
pard.roistartt.Width=.5;
pard.roistart.object=struct('Style','edit','String','7');
pard.roistart.position=[6,4.5];
pard.roistart.Width=.5;

pard.syncParameters={{'cal_3Dfile','cal_3Dfile',{'String'}}};

pard.plugininfo.type='WorkflowModule'; 
pard.plugininfo.description='Single-emiiter fitter that considers all ther candidate localizations. More robust for high-density localization than a single-emitter fitter';
end