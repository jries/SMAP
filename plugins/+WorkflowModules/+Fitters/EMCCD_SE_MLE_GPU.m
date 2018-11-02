classdef EMCCD_SE_MLE_GPU<interfaces.WorkflowFitter
    properties
        fitpar
        fitfunction
    end
    methods
        function obj=EMCCD_SE_MLE_GPU(varargin)
            obj@interfaces.WorkflowFitter(varargin{:})
            obj.inputChannels=1; 
             obj.setInputChannels(1,'frame');
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function fitinit(obj)
            obj.fitpar=getfitpar(obj);
            % check if x,y, then initialize range etc
            obj.fitfunction = @obj.nofound;
            disp('checking cuda version')
            reporttext='GPU fit function did not run. Possibly the wrong CUDA version is installed.';
            img=zeros(7,'single');img(3,3)=1;
            a50=0;a75=0;
            try 
                a50=gaussmlev2_cuda50(img,1,10,1);
            end
            try 
                a75=GPUgaussMLEv2_CUDA75(img,1,10,1);
            end
            
            if any(a75~=0)
                obj.fitfunction=@GPUgaussMLEv2_CUDA75;
                reporttext='GPUgaussMLEv2_CUDA75';  
            elseif any(a50~=0)
                  obj.fitfunction=@gaussmlev2_cuda50;
                reporttext='gaussmlev2_cuda50';                
            end
            
            roisize=obj.getPar('loc_ROIsize');
            obj.numberInBlock=round(5500*100/roisize^2);
            
            disp(reporttext)
            if exist('err','var')&&exist('err2','var')
                err
                err2
            end
%             try
%                reset(gpuDevice);
%             catch
%             end

        end
        function nofound(obj,varargin)
            disp('fit function not working. Wrong Cuda version?')
        end
        
        function locs=fit(obj,imstack,stackinfo)
             s=size(imstack);
             if length(s)==2 
                 s(3)=1;
             end
             if s(3)==0
                 locs=[];
                 return
             end
           fitpar=obj.fitpar;
            EMexcess=fitpar.EMexcessNoise;
            if isempty(EMexcess)
                EMexcess=1;
            end
           try %%for test set to 0: no fitting
               if obj.fitpar.fitmode==3
                   X=stackinfo.X;Y=stackinfo.Y;
                    zpar=[fitpar.zpar{X,Y}(:)];
%                     (ii{1},data,PSFSigma,iterations,fittype,Ax,Ay,Bx,By,gamma,d,PSFy0);
%                     [P CRLB LogL]=fitter.gaussmlev2_cuda70_newz(imstack,zpar(1),p.iterations,p.fitmode,zpar(2),zpar(3),zpar(4),zpar(5),zpar(6),zpar(7),zpar(8));
%                     [P CRLB LogL]=fitter.GPU_MLE_Lidke_55(imstack,zpar(1),p.iterations,p.fitmode,zpar(2),zpar(3),zpar(4)*0,zpar(5)*0,zpar(6),zpar(7),zpar(8));
%                      [P CRLB LogL]=fitter.GPUMLESINGEMITTER(imstack,zpar(1),p.iterations,p.fitmode,zpar(2),zpar(3),zpar(4)*0,zpar(5)*0,zpar(6),zpar(7),zpar(8));
                    [P CRLB LogL]=obj.fitfunction(imstack/EMexcess,zpar(1),obj.fitpar.iterations,obj.fitpar.fitmode,zpar(2),zpar(3),zpar(4),zpar(5),zpar(6),zpar(7),zpar(8));
               else
                   [P CRLB LogL]=obj.fitfunction(imstack/EMexcess,obj.fitpar.PSFx0,obj.fitpar.iterations,obj.fitpar.fitmode); 
%                     [P CRLB LogL]=fitter.gaussmlev2_cuda70_newz(imstack,p.PSFx0,p.iterations,p.fitmode);  
%                      [P CRLB LogL]=fitter.GPU_MLE_Lidke_55(imstack,p.PSFx0,p.iterations,p.fitmode);  
               end
           catch
               disp('gaussmlev2_cuda50 did not work')
               P=zeros(s(3),15);LogL=zeros(s(3),1);CRLB=P;
           end

           
           v1=ones(s(3),1,'single');
           
           dn=ceil((s(1)-1)/2)*v1;
           
            shiftx=0;%-0.5; %deviation from ground truth
             shifty=0;%-0.5;
             posx=stackinfo.x+shiftx;
             posy=stackinfo.y+shifty;
%            posx=stackinfo.x;
%            posy=stackinfo.y;
           frame=stackinfo.frame;
           locs.xpix=P(:,2)-dn+posx;
           locs.ypix=P(:,1)-dn+posy;
           locs.phot=P(:,3)*EMexcess;
           locs.bg=P(:,4)*EMexcess;
           locs.frame=frame;
           
           locs.xerrpix=sqrt(CRLB(:,2));
           locs.yerrpix=sqrt(CRLB(:,1));
           locs.photerr=sqrt(CRLB(:,3))*EMexcess;
           locs.bgerr=sqrt(CRLB(:,4))*EMexcess;
           locs.logLikelihood=LogL;
           
           locs.peakfindx=posx;
           locs.peakfindy=posy;

            switch obj.fitpar.fitmode
                case 1 %sx not fitted
                    sx=fitpar.PSFx0*v1;
                    locs.PSFxpix=0*locs.xpix+sx;
                    locs.PSFypix=locs.PSFxpix;
                case 2 % sx free
                    locs.PSFxpix=P(:,5);
                    locs.PSFxerr=sqrt(CRLB(:,5));
%                     sx=locs.PSFx;
                    locs.PSFypix=locs.PSFxpix;
                case 3
                    locs.znm=(P(:,5)*1000+obj.fitpar.objPos*v1)*fitpar.refractive_index_mismatch;
                    locs.zerr=sqrt(CRLB(:,5))*1000*fitpar.refractive_index_mismatch;
                    [locs.PSFxpix,locs.PSFypix]=zpar2sigma(locs.znm/1000,zpar);
                    

                case 4  %sx,sy
                    
                    locs.PSFxpix=P(:,5);
                    locs.PSFxerr=sqrt(CRLB(:,5));
                    locs.PSFypix=P(:,6);
                    locs.PSFyerr=sqrt(CRLB(:,6));  
%                     sx=locs.PSFx;
            end
            locs.locpthompson=sqrt((locs.PSFxpix.*locs.PSFypix+1/12*v1)./( locs.phot/EMexcess)+8*pi*(locs.PSFxpix.*locs.PSFypix).^2.* locs.bg./( locs.phot/EMexcess).^2);
            
        end

        
        function initGui(obj)
            initGui@interfaces.WorkflowFitter(obj);
            obj.guihandles.fitmode.Callback={@fitmode_callback,obj};
            fitmode_callback(0,0,obj)
            obj.guihandles.loadcal.Callback={@loadcall_callback,obj};

        end
%         function prerun(obj,p)
%             obj.fitpar=getfitpar(obj);   
%             prerun@interfaces.WorkflowFitter(obj);
%                   
%         end
        
            
    end
end

function loadcall_callback(a,b,obj)
p=obj.getAllParameters;
[f,p]=uigetfile(p.cal_3Dfile);
if f
    l=load([p f]);
    if ~isfield(l,'outforfit') && ~isfield(l,'SXY')
        msgbox('no 3D data recognized. Select other file.');
    end
    obj.setGuiParameters(struct('cal_3Dfile',[p f]));
    
end
end

function fitpar=getfitpar(obj)
p=obj.getAllParameters;
fitpar.iterations=p.iterations;
fitpar.fitmode=p.fitmode.Value;
if fitpar.fitmode==3
    calfile=p.cal_3Dfile;
    cal=load(calfile);
    if 0% p.useObjPos
        
        disp('obj. position not implemented yet')
    else
        fitpar.objPos=p.objPos;
        if isfield(cal,'outforfit')
            fitpar.zpar{1,1}=cal.outforfit;
        else
            s=size(cal.SXY);
            Z=1;
            if p.useObjPos
                zr=cal.SXY(1).Zrangeall;
                zr(1)=[];zr(end)=inf;
                Z=find(p.objPos<=zr,1,'first');
            end
            for X=1:s(1)
                for Y=1:s(2)
                    zpar{X,Y}=cal.SXY(X,Y,Z).fitzpar;
                end
            end
            fitpar.zpar=zpar;
            if size(cal.SXY,3)>1
                obj.spatial3Dcal=true;
            else
                obj.spatial3Dcal=false;
            end
            xr=cal.SXY(1,1).Xrangeall;
            xr(1)=-inf;xr(end)=inf;
            yr=cal.SXY(1,1).Yrangeall;
            yr(1)=-inf;yr(end)=inf;
            obj.spatialXrange=xr;
            obj.spatialYrange=yr;
                
        end
        fitpar.refractive_index_mismatch=p.refractive_index_mismatch;
    end
    
else
    fitpar.PSFx0=p.PSFx0;
end
if p.loc_cameraSettings.EMon
    fitpar.EMexcessNoise=2;
else
fitpar.EMexcessNoise=1;
end
% fitpar.EMexcessNoise
end

function fitmode_callback(a,b,obj)
p=obj.getGuiParameters;
fitmode=p.fitmode.Value;
switch fitmode
    case 3
        ton={'loadcal','cal_3Dfile','useObjPos','objPos','trefractive_index_mismatch','refractive_index_mismatch'};
        toff={'PSFx0','tPSFx0'};
    otherwise
        toff={'loadcal','cal_3Dfile','useObjPos','objPos','trefractive_index_mismatch','refractive_index_mismatch'};
        ton={'PSFx0','tPSFx0'};
end

switch fitmode
    case {1,2}
        roisize=7;
        iterations=50;
      
    otherwise
        roisize=11;
        iterations=200;
end

obj.setPar('loc_ROIsize',roisize);

obj.fieldvisibility('on',ton,'off',toff);
obj.setGuiParameters(struct('iterations',iterations));
% for k=1:length(ton)
%     
%     obj.guihandles.(ton{k}).Visible='on';
% end
% for k=1:length(toff)
%     obj.guihandles.(toff{k}).Visible='off';
% end
end

function pard=guidef
pard.fitmode.object=struct('Style','popupmenu','String','PSF fix|PSF free|3D z|ellipt: PSFx PSFy','Value',2);
pard.fitmode.position=[1,1];
pard.fitmode.Width=2;
pard.fitmode.TooltipString=sprintf('Fit mode. Fit with constant PSF, free PSF, 3D with astigmatism, asymmetric PSF (for calibrating astigmatic 3D)');

pard.text.object=struct('Style','text','String','Iterations:');
pard.text.position=[1,3.3];
pard.text.Width=0.7;
pard.text.Optional=true;
pard.iterations.object=struct('Style','edit','String','50');
pard.iterations.position=[1,4];
pard.iterations.TooltipString=sprintf('number of iterations for the GPU fitter (typical 50, use 100-250 for ellipt: PSFx PSFy or 3Dz).');
pard.iterations.Optional=true;

pard.tPSFx0.object=struct('Style','text','String','PSFx start (pix)');
pard.tPSFx0.position=[2,1];
pard.tPSFx0.Width=1.25;
pard.tPSFx0.Optional=true;

pard.PSFx0.object=struct('Style','edit','String','1');
pard.PSFx0.position=[2,2.25];
pard.PSFx0.Width=0.75;
pard.PSFx0.TooltipString=sprintf('start value for PSF, or size of PSF when PSF fixed (in camera pixels)');
pard.PSFx0.Optional=true;

pard.loadcal.object=struct('Style','pushbutton','String','Load 3D cal');
pard.loadcal.position=[2,1];
pard.cal_3Dfile.object=struct('Style','edit','String','settings/cal_3DAcal.mat');
pard.cal_3Dfile.position=[2,2];
pard.cal_3Dfile.Width=3;
pard.cal_3Dfile.TooltipString=sprintf('3D calibration file for astigmtic 3D. \n Generate from bead stacks with plugin: Analyze/sr3D/CalibrateAstig');

pard.useObjPos.object=struct('Style','checkbox','String','Use objective position:');
pard.useObjPos.position=[3,1];
pard.useObjPos.Width=3;
pard.useObjPos.Optional=true;

pard.objPos.object=struct('Style','edit','String','0');
pard.objPos.position=[3,4];
pard.objPos.TooltipString=sprintf('Position of the objective above the coverslip (nm, piezo position). \n Only used in combination with CalibrateAstigDeep.');
pard.objPos.Optional=true;

pard.trefractive_index_mismatch.object=struct('Style','text','String','Refractive Index mismatch factor:');
pard.trefractive_index_mismatch.position=[4,1];
pard.trefractive_index_mismatch.Width=3;
pard.trefractive_index_mismatch.Optional=true;

pard.refractive_index_mismatch.object=struct('Style','edit','String','.72');
pard.refractive_index_mismatch.position=[4,4];
pard.refractive_index_mismatch.TooltipString=sprintf('Correction factor to take into account the different refracrive indices of immersion oil and buffer. \n This leads to smaller distances inside the sample compared to bead calibration. \n Bead calibration: in piezo positions (nm). \n This factor transforms z positions to real-space z positions. \n For high-NA oil objectives: typical 0.72 (range 0.7-1).');
pard.refractive_index_mismatch.Optional=true;

pard.plugininfo.type='WorkflowFitter';
pard.plugininfo.description='Maximum likelyhood estimater, optimized for GPU processing. According to: C. S. Smith, N. Joseph, B. Rieger, and K. A. Lidke, ?Fast, single-molecule localization that achieves theoretically minimum uncertainty.,? Nat Methods, vol. 7, no. 5, pp. 373?375, May 2010.';
end