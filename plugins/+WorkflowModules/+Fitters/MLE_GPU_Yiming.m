classdef MLE_GPU_Yiming<interfaces.WorkflowFitter
    properties
        fitpar
%         mirrorstack
    end
    methods
        function obj=MLE_GPU_Yiming(varargin)
            obj@interfaces.WorkflowFitter(varargin{:})
            obj.inputChannels=1; 
             obj.setInputChannels(1,'frame');
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function fitinit(obj)
            obj.fitpar=getfitpar(obj);
            % check if x,y, then initialize range etc
            obj.fitpar.fitfunction = @obj.nofound;
%             disp('checking cuda fit')
%             reporttext='GPU fit function did not run. Possibly the wrong CUDA version is installed.';
            img=zeros(11,'single');img(5,5)=1;
            
%             try
%                 fitp=callYimingFitter(img,single(1),single(10),single(2),single(0),0);
                fitp=mleFit_LM(img,1);
                obj.fitpar.fitfunction=@mleFit_LM;
%                 obj.fitpar.fitfunction=@callYimingFitter;
%                 reporttext='mleFit_LM works';
%             end
            roisize=obj.getPar('loc_ROIsize');
            obj.numberInBlock=round(obj.fitpar.roisperfit*100/roisize^2);
            p=obj.getAllParameters;
            if obj.fitpar.fitmode==5 ||obj.fitpar.fitmode==6
%                 EMfile=obj.getPar('loc_fileinfo').EMon;
%                 EMcal=obj.fitpar.EMon;
%                 EMcal=obj.fitpar.splinefit{1}.isEM;
%                 mirrorstack=obj.getSingleGuiParameter('automirror');
%                 switch mirrorstack.selection
%                     case 'auto'
%                         obj.fitpar.mirrorstack=~(EMfile==EMcal);
%                     case 'mirror'
%                         obj.fitpar.mirrorstack=true;
%                     otherwise
%                         obj.fitpar.mirrorstack=false;
%                 end
%                 if obj.getSingleGuiParameter('automirror')
%                     
%                 else
%                     obj.fitpar.mirrorstack=false;
%                 end
                obj.fitpar.mirrorstack=false; %later: take out completely. NOw in tiff loader
                if p.overwritePixelsize
                    obj.setPar('overwrite_pixelsize',[p.pixelsizex p.pixelsizey])
                    cs=obj.getPar('loc_cameraSettings');
                    cs.cam_pixelsize_um=[p.pixelsizex p.pixelsizey];
                    obj.setPar('loc_cameraSettings',cs);
                else
                    obj.setPar('overwrite_pixelsize',[])
                end
                
            end   
            obj.setPar('loc_iterations',p.iterations);
%             disp(reporttext)
        end
        function nofound(obj,varargin)
            disp('fit function not working. Wrong Cuda version?')
        end

        function locs=fit(obj,imstack,stackinfo)
            if obj.fitpar.fitmode==3
                X=stackinfo.X;Y=stackinfo.Y;
                obj.fitpar.zparhere=[obj.fitpar.zpar{X,Y}(:)];
            elseif obj.fitpar.fitmode==5 || obj.fitpar.fitmode==6
                X=stackinfo.X;Y=stackinfo.Y;
                obj.fitpar.splinefithere=[obj.fitpar.splinefit{X,Y}(:)];
            end
            if obj.fitpar.issCMOS
                varstack=getvarmap(obj.fitpar.varmap,stackinfo,size(imstack,1));
%                 imstackraw=imstack; %XXXXXXXXXXXXXXXXXXXX
                if ~isempty(obj.fitpar.offsetmap)
                    offsetstack=getvarmap(obj.fitpar.offsetmap,stackinfo,size(imstack,1));
                    imstack=imstack-offsetstack;
                end
                if ~isempty(obj.fitpar.gainmap)
                    gainstack=getvarmap(obj.fitpar.gainmap,stackinfo,size(imstack,1));
                    imstack=imstack.*gainstack;
                end               
            else
                varstack=0;
            end
            out=fitwrapper(imstack,obj.fitpar,varstack);
            locs=fit2locs(out,stackinfo,obj.fitpar,imstack);
            
            if obj.fitpar.asymmetry
            [locs.asymmetry,locs.asymmdiag,locs.asymangle]=asymmetry(imgstack,true);
            end
            
            if ~isempty(locs)
                
            fn=fieldnames(stackinfo);
            infonames=setdiff(setdiff(fn,fieldnames(locs)), {'x','y','frame','Y','X'});
            locs=copyfields(locs,stackinfo,infonames);
            end
        end

        function initGui(obj)
            initGui@interfaces.WorkflowFitter(obj);
%             obj.guihandles.fitmode.Callback={@fitmode_callback,obj};
%             fitmode_callback(0,0,obj)
            obj.guihandles.loadcal.Callback={@loadcall_callback,obj};
%             obj.addSynchronization('loc_fileinfo',[],[],{@loc_fileinfo_callback,obj});

        end
            
    end
end

function loadscmos_callback(a,b,obj)
fs=obj.getSingleGuiParameter('scmosfile');
if isempty(fs)
    fs='*.*';
end
[file,pfad]=uigetfile(fs);
if file
    obj.setGuiParameters(struct('scmosfile',[pfad file]))
end

end

% function loc_fileinfo_callback(obj)
% 
% end
function locs=fit2locs(results,stackinfo,fitpar,image)
if isempty(results)
    locs=[];
    return
end
numl=size(results.P,1);

v1=ones(numl,1,'single');
s=size(image);          
dn=ceil((s(1)-1)/2)*v1;

shiftx=0;%-0.5; %deviation from ground truth
shifty=0;%-0.5;
posx=stackinfo.x+shiftx;
posy=stackinfo.y+shifty;
frame=stackinfo.frame;
P=results.P;
EMexcess=fitpar.EMexcessNoise;
CRLB=results.CRLB;
LogL=results.LogL;
           CRLB(isnan(CRLB))= 0; %XXXXXXXXX
           LogL(isnan(LogL))= 0; %XXXXXXXXX
           CRLB((CRLB)<0)= 0; %XXXXXXXXX
if (fitpar.fitmode==5||fitpar.fitmode==6) && fitpar.mirrorstack
    locs.xpix=dn-P(:,2)+posx;
else
    locs.xpix=P(:,2)-dn+posx;
end
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

switch fitpar.fitmode
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
        locs.znm=(P(:,5)*1000+fitpar.objPos*v1)*fitpar.refractive_index_mismatch;
        locs.zerr=sqrt(CRLB(:,5))*1000*fitpar.refractive_index_mismatch;
        [locs.PSFxpix,locs.PSFypix]=zpar2sigma(locs.znm/1000,fitpar.zparhere);


    case 4  %sx,sy

        locs.PSFxpix=P(:,5);
        locs.PSFxerr=sqrt(CRLB(:,5));
        locs.PSFypix=P(:,6);
        locs.PSFyerr=sqrt(CRLB(:,6));  
    case {5,6}
        
        %         locs.znm=(P(:,5)*1000+fitpar.objPos*v1)*fitpar.refractive_index_mismatch;
        locs.znm=-((P(:,5)-fitpar.z0)*fitpar.dz)*fitpar.refractive_index_mismatch;
%         notconverged=P(:,5)<2|P(:,5)>fitpar.coeffsize(3)-2;
%         locs.znm(notconverged)=NaN;
        
        locs.zerr=sqrt(CRLB(:,5))*fitpar.dz*fitpar.refractive_index_mismatch;
%         [locs.PSFxpix,locs.PSFypix]=zpar2sigma(locs.znm/1000,fitpar.zparhere);
        
         sx=1*v1;
        locs.PSFxpix=sx;
        locs.PSFypix=sx;
end
if fitpar.addgaussfit
    if fitpar.mirrorstack %now taken care of during loading
        xpixGauss=dn-results.P2(:,2)+posx;
    else
        xpixGauss=results.P2(:,2)-dn+posx;
    end
     ypixGauss=results.P2(:,1)-dn+posy;
     switch fitpar.addgaussfit_xyfrom.Value
         case 1 %spline
             locs.xpixGauss=xpixGauss;
             locs.ypixGauss=ypixGauss;
         case 2 %Gauss
             locs.xpixspline=locs.xpix;
             locs.ypixspline=locs.ypix;
             locs.xpix=xpixGauss;
             locs.ypix=ypixGauss;
             locs.xerrpix=sqrt(results.CRLB2(:,2));
             locs.yerrpix=sqrt(results.CRLB2(:,1)); 
     end
     switch fitpar.addgaussfit_mode.selection
         case 'free'
             locs.PSFxpixGauss=results.P2(:,5);
         case 'elliptical'
             locs.PSFxpixGauss=results.P2(:,5);
             locs.PSFypixGauss=results.P2(:,6);
     end
         
         
end
locs.locpthompson=sqrt((locs.PSFxpix.*locs.PSFypix+1/12*v1)./( locs.phot/EMexcess)+8*pi*(locs.PSFxpix.*locs.PSFypix).^2.* locs.bg./( locs.phot/EMexcess).^2);
locs.iterations=results.P(:,end);
end

function out=fitwrapper(imstack,fitpar,varstack)
s=size(imstack);
if length(s)==2 
 s(3)=1;
end
if s(3)==0
    out=[];
 return
end
% fitpar=obj.fitpar;
EMexcess=fitpar.EMexcessNoise;
if isempty(EMexcess)
    EMexcess=1;
end

arguments{2}=fitpar.fitmode;
arguments{3}=fitpar.iterations;

arguments{5}=varstack;
arguments{6}=1;
if isfield(fitpar,'dz')
    arguments{7}=fitpar.zstart/fitpar.dz;
end
    switch fitpar.fitmode
        case {1,2,4} %fix
            arguments{4}=fitpar.PSFx0;
            arguments{1}=single(imstack/EMexcess);
%         case 2 %free
        case 3 %z
            arguments{1}=single(imstack/EMexcess);
            arguments{4}=single(fitpar.zparhere);
%         case 4 %sx sy
        case {5,6} %spline   
            if fitpar.mirrorstack
                arguments{1}=single(imstack(:,end:-1:1,:)/EMexcess);
            else
                arguments{1}=single(imstack/EMexcess);
            end
            coeffh=(fitpar.splinefithere.cspline.coeff);
            if iscell(coeffh)
                coeffh=coeffh{1};
            end
            arguments{4}=single(coeffh);
            
    end
   
%     if fitpar.fitmode==6
%           
%         [P1 CRLB1 LL1 P2 CRLB2 LL2 ]=fitpar.fitfunction(arguments{:});
%         P = repmat(single(LL1>=LL2),1,6).*P1+repmat(single(LL1<LL2),1,6).*P2;
%         CRLB = repmat(single(LL1>=LL2),1,5).*CRLB1+repmat(single(LL1<LL2),1,5).*CRLB2;
%         LogL = repmat(single(LL1>=LL2),1,1).*LL1+repmat(single(LL1<LL2),1,1).*LL2;
%     else


        [P CRLB LogL]=fitpar.fitfunction(arguments{:});
        if fitpar.addgaussfit && fitpar.fitmode==5
            val2mode=[1 2 4];
            arguments{2}=val2mode(fitpar.addgaussfit_mode.Value);
            arguments{4}=1;
            [P2 CRLB2 LogL2]=fitpar.fitfunction(arguments{:});
        else
            P2=[];CRLB2=[];
        end
%     end

out.P=P;
out.CRLB=CRLB;
out.LogL=LogL;
out.P2=P2;
out.CRLB2=CRLB2;
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

function fitpar=getfitpar(obj)
p=obj.getAllParameters;
fitpar=p;
fitpar.addgaussfit=fitpar.addgaussfit & (fitpar.fitmode.Value ==5);
% fitpar.iterations=p.iterations;
fitpar.fitmode=p.fitmode.Value;
% fitpar.roisperfit=p.roisperfit;
fitpar.issCMOS=p.isscmos;
% fitpar.asymmetry=p.asymmetry;

if fitpar.fitmode==3||fitpar.fitmode==5
    fitpar.zstart=p.zstart;
  
    calfile=p.cal_3Dfile;
    cal=load(calfile);

    fitpar.objPos=0;
    if isfield(cal,'outforfit')
        fitpar.zpar{1,1}=cal.outforfit;
    elseif isfield(cal,'SXY')
        s=size(cal.SXY);
        Z=1;

        for X=s(1):-1:1
            for Y=s(2):-1:1
                zpar{X,Y}=cal.SXY(X,Y,Z).gauss_zfit;
                splinefit{X,Y}=cal.SXY(X,Y,Z);
            end
        end
        if ~isempty(splinefit{1})
            fitpar.dz=splinefit{1}.cspline.dz;
            fitpar.z0=splinefit{1}.cspline.z0;
            fitpar.splinefit=splinefit;
            coeffh=splinefit{1}.cspline.coeff;
            if iscell(coeffh)
                coeffh=coeffh{1};
            end
            fitpar.coeffsize=size(coeffh);
        end
        fitpar.zpar=zpar;

        if numel(cal.SXY)>1
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
        fitpar.EMon=cal.SXY(1).EMon;
    elseif isfield(cal,'cspline')
        fitpar.zpar{1}=cal.gauss_zfit;
        fitpar.dz=cal.cspline.dz;
        fitpar.z0=cal.cspline.z0;
        fitpar.splinefit{1}=cal.cspline_all;
        if ~isfield(fitpar.splinefit{1}.cspline,'isEM')
            fitpar.splinefit{1}.cspline.isEM=false;
        end
        fitpar.EMon=false;

    else
        disp('no calibration found')

    end
    
    if p.userefractive_index_mismatch
        fitpar.refractive_index_mismatch=p.refractive_index_mismatch;
    else
        fitpar.refractive_index_mismatch=1;
    end
    else
    fitpar.PSFx0=p.PSFx0;
end

%load sCMOS
if p.isscmos  %this needs to be extended. Include offset correction as well!
    [~,~,ext]=fileparts(p.scmosfile);
    switch ext
        case '.tif'
            varmaph=imread(p.scmosfile);
            offsetmap=[];
        case '.mat'
            l=load(p.scmosfile);
%             if isstruct(l)

                %fn=fieldnames(varmaph);
                %varmaph=varmaph.(fn{1});
                if isfield(l,'metadata')
                    metadata=l.metadata;
                else
                    metadata=p.loc_cameraSettings;
                end
                if isfield(l,'mean')
                offsetmaph=(single(l.mean)-metadata.offset)*metadata.pix2phot;
                else
                    offsetmaph=[];
                end
                if isfield(l,'gainmap')
                    gainmap=(single(l.gainmap))/metadata.pix2phot;
                else
                    gainmap=[];
                end
                varmaph=single(l.variance)*metadata.pix2phot^2;
%             end
        otherwise
            disp('could not load variance map. No sCMOS noise model used.')
            p.isscmos=false;
            fitpar.issCMOS=false;
            varmaph=[];
            offsetmaph=[];
    end
    if ~isempty(varmaph)
        roi=p.loc_cameraSettings.roi;
%         varmap=varmaph(max(1,roi(1)):roi(3),max(1,roi(2)):roi(4)); %or +1? roi: width height or coordinates? This is wrong
        varmap=varmaph(roi(1)+1:roi(1)+roi(3),roi(2)+1:roi(2)+roi(4)); 
    end
    if ~isempty(offsetmaph)
        roi=p.loc_cameraSettings.roi;
%         varmap=varmaph(max(1,roi(1)):roi(3),max(1,roi(2)):roi(4)); %or +1? roi: width height or coordinates? This is wrong
        offsetmap=offsetmaph(roi(1)+1:roi(1)+roi(3),roi(2)+1:roi(2)+roi(4)); 
    end
    if ~isempty(gainmap)
        roi=p.loc_cameraSettings.roi;
%         varmap=varmaph(max(1,roi(1)):roi(3),max(1,roi(2)):roi(4)); %or +1? roi: width height or coordinates? This is wrong
        gainmap=gainmap(roi(1)+1:roi(1)+roi(3),roi(2)+1:roi(2)+roi(4)); 
    end
else 
    varmap=[];
    offsetmap=[];
    gainmap=[];
end
% fitpar.varmap=varmap*p.loc_cameraSettings.pix2phot;
fitpar.varmap=varmap;
fitpar.offsetmap=offsetmap; % added by Robin
fitpar.gainmap=gainmap; % added by Robin



if p.loc_cameraSettings.EMon
    fitpar.EMexcessNoise=2;
else
    fitpar.EMexcessNoise=1;
end

end

function varstack=getvarmap(varmap,stackinfo,roisize)
numim=length(stackinfo.x);
varstack=zeros(roisize,roisize,numim,'single');
dn=floor(roisize/2);
for k=1:numim
%     stackinfo.x(k)
%     stackinfo.y(k)
  varstack(:,:,k)=varmap(stackinfo.y(k)-dn:stackinfo.y(k)+dn,stackinfo.x(k)-dn:stackinfo.x(k)+dn);
    
end
end

function fitmode_callback(a,b,obj)
p=obj.getGuiParameters;
fitmode=p.fitmode.Value;
% fitz={'loadcal','cal_3Dfile','trefractive_index_mismatch','refractive_index_mismatch','overwritePixelsize','pixelsizex','pixelsizey','automirror','fit2D','isscmos','selectscmos','scmosfile'};
% fitxy={'PSFx0','tPSFx0'};
% switch fitmode
%     case {3,5}
%         ton=fitz;
%         toff=fitxy;
%     otherwise
%         toff=fitz;
%         ton=fitxy;
% end

switch fitmode
    case {1,2}
        roisize=7;
        iterations=50;
      
    otherwise
        roisize=13;
        iterations=50;
end

obj.setPar('loc_ROIsize',roisize);

% obj.fieldvisibility('on',ton,'off',toff);
obj.setGuiParameters(struct('iterations',iterations));
end

function pard=guidef(obj)
p1(1).value=1; p1(1).on={'PSFx0','tPSFx0'}; 
p1(1).off={'loadcal','cal_3Dfile','trefractive_index_mismatch','refractive_index_mismatch','overwritePixelsize',...
    'automirror','fit2D','pixelsizex','pixelsizey','addgaussfit','addgaussfit_mode','addgaussfit_t',...
    'addgaussfit_xyfrom'};%,'isscmos','selectscmos','scmosfile'};
p1(2)=p1(1);p1(2).value=2;
p1(3).value=3;p1(3).off={'PSFx0','tPSFx0'};p1(3).on={'loadcal','cal_3Dfile',...
    'trefractive_index_mismatch','refractive_index_mismatch',...
    'overwritePixelsize','automirror','fit2D'};%,'isscmos'};
p1(4)=p1(1);p1(4).value=4;
p1(5)=p1(3);p1(5).value=5;
p1(6)=p1(5);p1(6).value=6;
p1(5).on=[p1(5).on {'addgaussfit'}];



pard.fitmode.object=struct('Style','popupmenu','String',{{'PSF fix','PSF free','3D z','ellipt: PSFx PSFy','Spline'}},'Value',2,'Callback',{{@obj.switchvisible,p1,{@fitmode_callback,0,0,obj}}});
pard.fitmode.position=[1,1];
pard.fitmode.Width=1.5;
pard.fitmode.TooltipString=sprintf('Fit mode. Fit with constant PSF, free PSF, 3D with astigmatism, asymmetric PSF (for calibrating astigmatic 3D)');

pard.text.object=struct('Style','text','String','Iterations:');
pard.text.position=[1,2.5];
pard.text.Width=0.7;
pard.text.Optional=true;
pard.iterations.object=struct('Style','edit','String','50');
pard.iterations.position=[1,3.2];
pard.iterations.TooltipString=sprintf('number of iterations for the GPU fitter (typical 50, use 100-250 for ellipt: PSFx PSFy or 3Dz).');
pard.iterations.Optional=true;
pard.iterations.Width=0.5;

pard.roisperfitt.object=struct('Style','text','String','ROIs/fit:');
pard.roisperfitt.position=[1,3.9];
pard.roisperfitt.Width=0.6;
pard.roisperfitt.Optional=true;
pard.roisperfit.object=struct('Style','edit','String','15000');
pard.roisperfit.position=[1,4.5];
pard.roisperfit.TooltipString=sprintf('Number of 10 x 10 pixel ROIs passed to GPU for fitting. For other ROI sizes, the number is adjusted accordingly.');
pard.roisperfit.Optional=true;
pard.roisperfitt.TooltipString=pard.roisperfit.TooltipString;
pard.roisperfit.Width=0.5;


p(1).value=0; p(1).on={}; p(1).off={'addgaussfit_mode','addgaussfit_t','addgaussfit_xyfrom'};
p(2).value=1; p(2).on=p(1).off; p(2).off={};
pard.addgaussfit.object=struct('Style','checkbox','String','additional Gauss fit','Callback',{{@obj.switchvisible,p}});
pard.addgaussfit.position=[3,1];
pard.addgaussfit.Optional=true;
pard.addgaussfit.Width=1.3;
pard.addgaussfit_mode.object=struct('Style','popupmenu','String',{{'fix','free','elliptical'}},'Value',2);
pard.addgaussfit_mode.position=[3,2.3];
pard.addgaussfit_mode.Width=0.7;
pard.addgaussfit_mode.Optional=true;
pard.addgaussfit_t.object=struct('Style','text','String','x,y coordinates from');
pard.addgaussfit_t.position=[3,3];
pard.addgaussfit_t.Width=1.3;
pard.addgaussfit_t.Optional=true;
pard.addgaussfit_xyfrom.object=struct('Style','popupmenu','String',{{'spline','Gauss'}});
pard.addgaussfit_xyfrom.position=[3,4.1];
pard.addgaussfit_xyfrom.Optional=true;
pard.addgaussfit_xyfrom.Width=.9;

pard.tPSFx0.object=struct('Style','text','String','PSFx start (pix)');
pard.tPSFx0.position=[2,1];
pard.tPSFx0.Width=1.25;
pard.tPSFx0.Optional=true;

pard.PSFx0.object=struct('Style','edit','String','1');
pard.PSFx0.position=[2,2.25];
pard.PSFx0.Width=0.75;
pard.PSFx0.TooltipString=sprintf('start value for PSF, or size of PSF when PSF fixed (in camera pixels)');
pard.PSFx0.Optional=true;

% pard.fit2D.object=struct('Style','checkbox','String','2D PSF','Value',0);
% pard.fit2D.position=[2,3.5];
% pard.fit2D.Width=.75;
% pard.fit2D.TooltipString=sprintf('Check if PSF model is 2D (no specific PSF engineering), or displays a high degree of similarity above and below the focal plane');
% pard.fit2D.Optional=true;

pard.zstartt.object=struct('Style','text','String','z start (nm)');
pard.zstartt.position=[2,3.5];
pard.zstartt.Width=.75;
pard.zstartt.TooltipString=sprintf('z start parameter. Use vector with several values for 2D PSF');
pard.zstartt.Optional=true;

pard.zstart.object=struct('Style','edit','String','0');
pard.zstart.position=[2,4.25];
pard.zstart.Width=.75;
pard.zstart.TooltipString=pard.zstartt.TooltipString;
pard.zstart.Optional=true;

% pard.automirror.object=struct('Style','popupmenu','String',{{'auto','no mirror','mirror'}});
% pard.automirror.position=[2,4.25];
% pard.automirror.Width=.75;
% pard.automirror.Optional=true;

pard.loadcal.object=struct('Style','pushbutton','String','Load 3D cal');
pard.loadcal.position=[2,1];
pard.loadcal.Width=.75;
pard.cal_3Dfile.object=struct('Style','edit','String','');
pard.cal_3Dfile.position=[2,1.75];
pard.cal_3Dfile.Width=1.75;
pard.cal_3Dfile.TooltipString=sprintf('3D calibration file for astigmtic 3D. \n Generate from bead stacks with plugin: Analyze/sr3D/CalibrateAstig');


% p(1).value=0; p(1).on={}; p(1).off={'linkt','link'};
% p(2).value=1; p(2).on={'linkt','link'}; p(2).off={};
% 
% pard.isglobal.object=struct('Style','checkbox','String','Global fit','Callback',{{@obj.switchvisible,p}});
% pard.isglobal.position=[3,1];
% pard.isglobal.Width=1;
% pard.isglobal.Optional=true;
% 
% pard.linkt.object=struct('Style','text','String','link: x y z N bg');
% pard.linkt.position=[3,2];
% pard.linkt.Width=1.5;
% pard.linkt.Optional=true;
% 
% pard.link.object=struct('Style','edit','String','1 1 1 0 0');
% pard.link.position=[3,3.5];
% pard.link.Width=1.5;
% pard.link.Optional=true;

p(1).value=0; p(1).on={}; p(1).off={'refractive_index_mismatch'};
p(2).value=1; p(2).on={'refractive_index_mismatch'}; p(2).off={};
pard.userefractive_index_mismatch.object=struct('Style','checkbox','String','RI mismatch:','Callback',{{@obj.switchvisible,p}});
pard.userefractive_index_mismatch.position=[4,3.5];
pard.userefractive_index_mismatch.Width=1.5;
pard.userefractive_index_mismatch.Optional=true;

pard.refractive_index_mismatch.object=struct('Style','edit','String','.8');
pard.refractive_index_mismatch.position=[4,4.5];
pard.refractive_index_mismatch.TooltipString=sprintf('Correction factor to take into account the different refracrive indices of immersion oil and buffer. \n This leads to smaller distances inside the sample compared to bead calibration. \n Bead calibration: in piezo positions (nm). \n This factor transforms z positions to real-space z positions. \n For high-NA oil objectives: typical 0.72 (range 0.7-1).');
pard.refractive_index_mismatch.Optional=true;
pard.refractive_index_mismatch.Width=0.5;


p(1).value=0; p(1).on={}; p(1).off={'pixelsizex','pixelsizey'};
p(2).value=1; p(2).on={'pixelsizex','pixelsizey'}; p(2).off={};
pard.overwritePixelsize.object=struct('Style','checkbox','String','New pixelsize X,Y (um):','Callback',{{@obj.switchvisible,p}});
pard.overwritePixelsize.position=[4,1];
pard.overwritePixelsize.Width=1.5;
pard.overwritePixelsize.Optional=true;

pard.pixelsizex.object=struct('Style','edit','String','.1');
pard.pixelsizex.position=[4,2.5];
pard.pixelsizex.Width=0.5;
pard.pixelsizex.Optional=true;

pard.pixelsizey.object=struct('Style','edit','String','.1');
pard.pixelsizey.position=[4,3];
pard.pixelsizey.Width=0.5;
pard.pixelsizey.Optional=true;

p(1).value=0; p(1).on={}; p(1).off={'selectscmos','scmosfile'};
p(2).value=1; p(2).on={'selectscmos','scmosfile'}; p(2).off={};
pard.isscmos.object=struct('Style','checkbox','String','sCMOS','Callback',{{@obj.switchvisible,p}});   
pard.isscmos.position=[5,1];
pard.isscmos.Optional=true;
pard.selectscmos.object=struct('Style','pushbutton','String','Load var map','Callback',{{@loadscmos_callback,obj}});   
pard.selectscmos.TooltipString='Select sCMOS variance map (in ADU^2) of same size ROI on chip as image stack';
pard.selectscmos.position=[5,2];
pard.selectscmos.Optional=true;
pard.scmosfile.object=struct('Style','edit','String','');
pard.scmosfile.TooltipString='Tiff/.mat image containing sCMOS variance map (same ROI on camera as tiff).';
pard.scmosfile.position=[5,3];
pard.scmosfile.Optional=true;
    pard.scmosfile.Width=2;
   
pard.asymmetry.object=struct('Style','checkbox','String','get asymmetry');   
pard.asymmetry.position=[6,1];
pard.asymmetry.Optional=true;

pard.plugininfo.type='WorkflowFitter';
pard.plugininfo.description='Maximum likelyhood estimater, optimized for GPU processing. According to: C. S. Smith, N. Joseph, B. Rieger, and K. A. Lidke, ?Fast, single-molecule localization that achieves theoretically minimum uncertainty.,? Nat Methods, vol. 7, no. 5, pp. 373?375, May 2010.';
end