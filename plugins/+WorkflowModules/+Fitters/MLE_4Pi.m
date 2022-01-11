classdef MLE_4Pi<interfaces.WorkflowFitter
%     Global fitter for multiple channels.
    properties
        fitpar
        varmap
    end
    methods
        function obj=MLE_4Pi(varargin)
            obj@interfaces.WorkflowFitter(varargin{:})
            obj.inputChannels=1; 
             obj.setInputChannels(1,'frame');
             
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function fitinit(obj)
            p=obj.getAllParameters;
            obj.infofields={'xpix','ypix','ID','dx','dy'};
            obj.fitpar=getfitpar(obj);
            obj.fitpar.link=obj.fitpar.link([2 1 4 5 3 6]);
%             obj.fitpar.fitfunction=@mleFit_LM_4Pi;
%             obj.fitpar.transform=obj.getPar('loc_globaltransform');
            roisize=obj.getPar('loc_ROIsize');
            obj.fitpar.Boxsize=roisize;
            obj.numberInBlock=round(obj.fitpar.roisperfit*100/roisize^2/12)*12;
            obj.setPar('loc_iterations',p.iterations);
            obj.setPar('loc_numberOfChannels',4);
        end
        function nofound(obj,varargin)
            disp('fit function not working. Wrong Cuda version?')
        end

        function locs=fit(obj,imstack,stackinfo)
            if obj.fitpar.isscmos
                varstack=getvarmap(obj,stackinfo,size(imstack,1));
            else
                varstack=0;
            end
            out=fitwrapper_4pi(imstack,obj.fitpar,stackinfo,varstack);
            locs=fit2locs_4pi(out,stackinfo,obj.fitpar,imstack);
        end

        function initGui(obj)
            initGui@interfaces.WorkflowFitter(obj);
            obj.guihandles.loadcal.Callback={@loadcall_callback,obj};

            %make global fit control
            hold=obj.guihandles.globaltable;
            gt=uitable(hold.Parent,'Position',hold.Position);
            
            gt.RowName={'x','y','z','N','Bg','ph'};
            gt.ColumnName={'l','x'};
            
            hh=obj.guiPar.fontsize+5;
            wh=hold.Position(3)-hh-50;
            gt.ColumnWidth={hh,wh};
            gt.Data={true,'1';true,'1';true,'1';false,'1';false,'1';true,'1'};
            gt.ColumnEditable=true;
            obj.guihandles.globaltable=gt;
            hold.Visible='off';
            hold.delete;

        end
            
    end
end


function locs=fit2locs_4pi(results,stackinfo,fitpar,image)
if isempty(results)
    locs=[];
    return
end
numl=size(results.P,1);
numpar=6;
numchannels=4;

v1=ones(numl,1,'single');
s=size(image);          
dn=ceil((s(1)-1)/2)*v1;

shiftx=0;%-0.5; %deviation from ground truth
shifty=0;%-0.5;
posx=stackinfo(1).xpix+shiftx;
posy=stackinfo(1).ypix+shifty;
frame=stackinfo(1).frame;
P=results.P;
EMexcess=1;
CRLB=results.CRLB;
LogL=results.LogL;
           CRLB(isnan(CRLB))= 0; %XXXXXXXXX
           LogL(isnan(LogL))= 0; %XXXXXXXXX
           CRLB((CRLB)<0)= 0; %XXXXXXXXX
          % LogL((LogL)<0)= 0; %XXXXXXXXX
           
o=ones( numpar,1);fac=num2cell(o);z=zeros( numpar,1);off=num2cell(z);
faccrlb=fac;

fac{2}=1;
off{2}=-dn+posx+1; %+1 added to overlay raw and image

off{1}=-dn+posy+1;

fac{4}=EMexcess;
fac{3}=EMexcess;
faccrlb{4}=EMexcess;
faccrlb{3}=EMexcess;
locs.frame=frame;

locs.logLikelihood=LogL;

locs.peakfindx=posx;
locs.peakfindy=posy;
sx=1*v1; 
locs.PSFxpix=sx;
locs.PSFypix=sx;
fitpar.refractive_index_mismatch=1;

fac{5}=-fitpar.Pixelsizez*fitpar.refractive_index_mismatch; %fix direction
off{5}=+fitpar.Zcenter*fitpar.Pixelsizez*fitpar.refractive_index_mismatch;
faccrlb{5}=fitpar.Pixelsizez*fitpar.refractive_index_mismatch;

names={'ypix','xpix','phot','bg','zastig','phase'};
namesav={'ypix','xpix','zastig','phot'};
linked=fitpar.link(1:6);

ind=1;
for k=1:length(names)
    if linked(k)
        locs.(names{k})=P(:,ind).*fac{k}+off{k};
        locs.([names{k} 'err'])=sqrt(CRLB(:,ind)).*faccrlb{k};
        ind=ind+1;
    else
        v=zeros(numl,numchannels,'single');
        ve=zeros(numl,numchannels,'single');
        for c=1:numchannels
            ch=num2str(c);
            locs.([names{k} ch])=P(:,ind).*fac{k}+off{k};
            if strcmp(names{k},'phase')
                locs.([names{k} ch])=mod(locs.([names{k} ch]),2*pi);
            end
            locs.([names{k} ch 'err'])=sqrt(CRLB(:,ind)).*faccrlb{k};
            v(:,c)=locs.([names{k} ch]);
            ve(:,c)=locs.([names{k} ch 'err']);
            ind=ind+1;
        end
        if any(contains(namesav,names{k}))
%             switch fitpar.mainchannel
%                 case 1
                    locs.(names{k})=sum(v./ve,2)./sum(1./ve,2);
                    locs.([names{k} 'err'])=1./sqrt(1./ve(:,1).^2+1./ve(:,2).^2);
%                 case 2
%                     locs.(names{k})=v(:,1);
%                     locs.([names{k} 'err'])=ve(:,1);
%                 case 3
%                     locs.(names{k})=v(:,2);
%                     locs.([names{k} 'err'])=ve(:,2);
%             end

        end
    end
end
locs.iterations=P(:,end);
end

function out=fitwrapper_4pi(imstack,fitpar,stackinfo,varstack)
% numberOfChannels=4;
% nfits=ceil(size(imstack,3));
% npar=6;
% numframes=size(imstack,3); 

% s=size(imstack);
dx=[stackinfo(:).dx];
dy=[stackinfo(:).dy];
[P,CRLB,LogL] = psfloc_4Pi(fitpar,imstack,dy,dx,fitpar.link);

% EMexcess=1;


% if isfield(fitpar.splinefithere.coeff,'normf')
%     normf=fitpar.splinefithere.coeff.normf;
% else
%     normf=ones(4,1);
% end
        
% dT=zeros(npar,numberOfChannels,(nfits),'single');
% % imfit=zeros(s(1),s(2),ceil(s(3)),numberOfChannels,'single');
% for k=1:numberOfChannels
%     dT(1,k,:)=stackinfo(k).dy;
%     dT(2,k,:)=stackinfo(k).dx;
% %     if fitpar.mirror{k}==1 || fitpar.mirror{k}==3  %mirror x or xy
% %         imfit(:,:,:,k)=imstack(end:-1:1,:,k:numberOfChannels:end)/normf(k);
% %         dT(1,k,:)=-dT(1,k,:);
% %     else
% %         imfit(:,:,:,k)=imstack(:,:,k:numberOfChannels:end)/normf(k);
% %     end
%     
% end
% sharedA = int32(fitpar.link);
% sharedA(end+1)=1; %link phase. Not part of GUI
% 
% out.indused=1:numberOfChannels:numframes;

% [P,CRLB1 LL] = CPUmleFit_LM_MultiChannel(single(imstack(:, :, :, :)),uint32(shared),iterations,single(PSF.Ispline), single(PSF.Aspline),single(PSF.Bspline),single(dTAll),single(phi0),z0);
% imageslicer(residual)
% 
% arguments{1}=single(imstack);
% arguments{2}=uint32(sharedA);
% arguments{3}=fitpar.iterations;
% arguments{4}=single(fitpar.splinefithere.coeff.Ispline);
% arguments{5}=single(fitpar.splinefithere.coeff.Aspline);
% arguments{6}=single(fitpar.splinefithere.coeff.Bspline);
% arguments{7}=single(dT); %XXXXXXXXXX
% %imstack, sharedflag, iterations, spline coefficients, channelshift,
% %fitmode, varmap
% arguments{8}=single(fitpar.splinefithere.coeff.phaseshifts);
% arguments{9}=fitpar.zstart/fitpar.splinefithere.dz+fitpar.splinefithere.z0;
% [P CRLB LogL]=fitpar.fitfunction(arguments{:});

out.P=P;
out.CRLB=CRLB;
out.LogL=LogL;
end
        

function loadcall_callback(a,b,obj)
p=obj.getAllParameters;
ph=fileparts(p.cal_3Dfile);
if ~exist(ph,'file')
    p.cal_3Dfile= [fileparts(obj.getPar('loc_outputfilename')) filesep '*.mat'];
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

function fitpar=getfitpar(obj)
%get all necessary parameters for fitting and store them 
p=obj.getAllParameters;

% fitpar.mainchannel=p.mainchannel.Value;
% fitpar.weightsch=p.weightsch;
% fitpar.fixPhot=p.fixPhot;
% if fitpar.fixPhot
%     fitpar.PhotonRatios=p.PhotonRatios;
% else
%     fitpar.PhotonRatios=[];
% end
% if length(fitpar.weightsch)==1
%     fitpar.weightsch(2)=1;
% end
% for k=1:size(p.globaltable.Data,1)
%     fitpar.factor{k}=str2num(p.globaltable.Data{k,2});
% end
calfile=p.cal_3Dfile;
cal=load(calfile);
fitpar=resetparam_4Pi_fit('',cal,0,[0 0]);
fitpar.iterations=p.iterations;
fitpar.roisperfit=p.roisperfit;
fitpar.isscmos=p.isscmos;
fitpar.link=[p.globaltable.Data{:,1}];
fitpar.zstart=p.zstart;
fitpar.isscmos= p.isscmos;
  
fitpar.mode='4pi';

% fitpar.dz=cal.pixelsize_z;
% fitpar.pixelsize=cal.pixelsize_x;
% fitpar.zT=cal.zT;

savefit=copyfields([],cal,{'T','zT','pixelsize_x','pixelsize_z','bead_radius'});
savefit.cal_3Dfile=p.cal_3Dfile;
obj.setPar('savefit',struct('cal3D',savefit));

% fitpar.splinecoeff=permute(cal.coeff,6:-1:1);

    %load sCMOS
%     if p.isscmos
%         varmaph=obj.getPar('cam_varmap');
%         if ~isempty(varmaph)
%             roi=p.loc_cameraSettings.roi;
%             varmap=varmaph(roi(1)+1:roi(1)+roi(3),roi(2)+1:roi(2)+roi(4)); 
%         else 
%             varmap=[];
%         end
%     else 
%         varmap=[];
%     end
%     fitpar.varmap=varmap*p.loc_cameraSettings.pix2phot;

% fitpar.EMexcessNoise=1;
end

function varstack=getvarmap(obj,stackinfo,roisize)
if isempty(obj.varmap)
    obj.varmap=obj.getPar('cam_varmap');
    if isempty(obj.varmap)
        disp('no sCMOS variance map found');
    end
end
numim=length(stackinfo.xpix);
varstack=zeros(roisize,roisize,numim,'single');
dn=floor(roisize/2);
for k=1:numim
    varstack(:,:,k)=obj.varmap(stackinfo.ypix(k)-dn:stackinfo.ypix(k)+dn,stackinfo.xpix(k)-dn:stackinfo.xpix(k)+dn);
    %varstack(:,:,k)=obj.varmap(stackinfo.xpix(k)-dn:stackinfo.xpix(k)+dn,stackinfo.ypix(k)-dn:stackinfo.ypix(k)+dn);
end
end


function pard=guidef(obj)
p1(1).value=1; p1(1).on={'PSFx0','tPSFx0'}; 
p1(1).off={'loadcal','cal_3Dfile','userefractive_index_mismatch','refractive_index_mismatch','overwritePixelsize',...
    'fit2D','pixelsizex','pixelsizey'};
p1(2).value=2;p1(2).off={'PSFx0','tPSFx0'};p1(2).on={'loadcal','cal_3Dfile','userefractive_index_mismatch','refractive_index_mismatch','overwritePixelsize','fit2D'};


% pard.fitmode.object=struct('Style','popupmenu','String',{{'PSF free','Spline'}},'Value',1,'Callback',{{@obj.switchvisible,p1,{@fitmode_callback,0,0,obj}}});
% pard.fitmode.position=[1,1];
% pard.fitmode.Width=1.5;
% pard.fitmode.TooltipString=sprintf('Fit mode. Fit with constant PSF, free PSF, 3D with astigmatism, asymmetric PSF (for calibrating astigmatic 3D)');

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
% 
% pard.tPSFx0.object=struct('Style','text','String','PSFx start (pix)');
% pard.tPSFx0.position=[2,1];
% pard.tPSFx0.Width=1.25;
% pard.tPSFx0.Optional=true;

% pard.PSFx0.object=struct('Style','edit','String','1');
% pard.PSFx0.position=[2,2.25];
% pard.PSFx0.Width=0.75;
% pard.PSFx0.TooltipString=sprintf('start value for PSF, or size of PSF when PSF fixed (in camera pixels)');
% pard.PSFx0.Optional=true;

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
pard.cal_3Dfile.object=struct('Style','edit','String','settings/cal_3DAcal.mat');
pard.cal_3Dfile.position=[2,1.75];
pard.cal_3Dfile.Width=1.75;
pard.cal_3Dfile.TooltipString=sprintf('3D calibration file for astigmtic 3D. \n Generate from bead stacks with plugin: Analyze/sr3D/CalibrateAstig');


% p(1).value=0; p(1).on={}; p(1).off={'globaltable','linkt','link','mainchannelt','mainchannel','weightch2t','weightch2'};
% p(2).value=1; p(2).on={'globaltable','linkt','link','mainchannelt','mainchannel','weightch2t','weightch2'}; p(2).off={};
% 
% pard.isglobal.object=struct('Style','checkbox','String','Global fit','Callback',{{@obj.switchvisible,p}});
% pard.isglobal.position=[3,3.];
% pard.isglobal.Width=.75;
% pard.isglobal.Optional=true;


% pard.mainchannelt.object=struct('Style','text','String','main x,y:');
% pard.mainchannelt.position=[4,3.];
% pard.mainchannelt.Width=.75;
% pard.mainchannelt.Optional=true;
% 
% pard.mainchannel.object=struct('Style','popupmenu','String',{{'mean','ch1','ch2'}});
% pard.mainchannel.position=[4,3.5];
% pard.mainchannel.Width=.75;
% pard.mainchannel.Optional=true;

% pard.weightscht.object=struct('Style','text','String','Weights r t');
% pard.weightscht.position=[5,3.];
% pard.weightscht.Width=.75;
% pard.weightscht.Optional=true;
% pard.weightsch.object=struct('Style','edit','String','1 1');
% pard.weightsch.position=[5,3.75];
% pard.weightsch.Width=.5;
% pard.weightsch.Optional=true;


pard.globaltable.object=struct('Style','listbox','String','x');
pard.globaltable.position=[7,4.25];
pard.globaltable.Width=.75;
pard.globaltable.Optional=true;
pard.globaltable.Height=5;
% pard.linkt.object=struct('Style','text','String','link: x y z N bg');
% pard.linkt.position=[3,2];
% pard.linkt.Width=1.5;
% pard.linkt.Optional=true;

% pard.link.object=struct('Style','edit','String','1 1 1 0 0');
% pard.link.position=[3,3.5];
% pard.link.Width=1.5;
% pard.link.Optional=true;
% p(1).value=0; p(1).on={}; p(1).off={'refractive_index_mismatch'};
% p(2).value=1; p(2).on={'refractive_index_mismatch'}; p(2).off={};
% pard.userefractive_index_mismatch.object=struct('Style','checkbox','String','RI mismatch:','Callback',{{@obj.switchvisible,p}});
% pard.userefractive_index_mismatch.position=[3,1];
% pard.userefractive_index_mismatch.Width=1.25;
% pard.userefractive_index_mismatch.Optional=true;


% pard.refractive_index_mismatch.object=struct('Style','edit','String','.8');
% pard.refractive_index_mismatch.position=[3,2.25];
% pard.refractive_index_mismatch.TooltipString=sprintf('Correction factor to take into account the different refracrive indices of immersion oil and buffer. \n This leads to smaller distances inside the sample compared to bead calibration. \n Bead calibration: in piezo positions (nm). \n This factor transforms z positions to real-space z positions. \n For high-NA oil objectives: typical 0.72 (range 0.7-1).');
% pard.refractive_index_mismatch.Optional=true;
% pard.refractive_index_mismatch.Width=0.35;

% 
% p(1).value=0; p(1).on={}; p(1).off={'pixelsizex','pixelsizey'};
% p(2).value=1; p(2).on={'pixelsizex','pixelsizey'}; p(2).off={};
% pard.overwritePixelsize.object=struct('Style','checkbox','String','pixelsize X,Y (um):','Callback',{{@obj.switchvisible,p}});
% pard.overwritePixelsize.position=[4,1];
% pard.overwritePixelsize.Width=1.25;
% pard.overwritePixelsize.Optional=true;

% pard.pixelsizex.object=struct('Style','edit','String','.1');
% pard.pixelsizex.position=[4,2.25];
% pard.pixelsizex.Width=0.35;
% pard.pixelsizex.Optional=true;
% 
% pard.pixelsizey.object=struct('Style','edit','String','.1');
% pard.pixelsizey.position=[4,2.6];
% pard.pixelsizey.Width=0.35;
% pard.pixelsizey.Optional=true;

% p(1).value=0; p(1).on={}; p(1).off={'selectscmos','scmosfile'};
% p(2).value=1; p(2).on={'selectscmos','scmosfile'}; p(2).off={};
pard.isscmos.object=struct('Style','checkbox','String','sCMOS');%,'Callback',{{@obj.switchvisible,p}});   
pard.isscmos.position=[7,1];
pard.isscmos.Optional=true;
% pard.selectscmos.object=struct('Style','pushbutton','String','Load var map','Callback',{{@loadscmos_callback,obj}});   
% pard.selectscmos.TooltipString='Select sCMOS variance map (in ADU^2) of same size ROI on chip as image stack';
% pard.selectscmos.position=[7,1.6];
% pard.selectscmos.Optional=true;
% pard.selectscmos.Width=.9;
% pard.scmosfile.object=struct('Style','edit','String','');
% pard.scmosfile.TooltipString='Tiff/.mat image containing sCMOS variance map (same ROI on camera as tiff).';
% pard.scmosfile.position=[7,2.5];
% pard.scmosfile.Optional=true;
%     pard.scmosfile.Width=.5;
    
% pard.asymmetry.object=struct('Style','checkbox','String','get asymmetry');   
% pard.asymmetry.position=[7,3];
% pard.asymmetry.Optional=true;
    
% 
% p(1).value=0; p(1).on={}; p(1).off={'PhotonRatios'};
% p(2).value=1; p(2).on={'PhotonRatios'}; p(2).off={};
% 
% pard.fixPhot.object=struct('Style','checkbox','String','Multi color: fix ratio to: ','Callback',{{@obj.switchvisible,p}});   
% pard.fixPhot.position=[6,1];
% pard.fixPhot.Optional=true;
% pard.fixPhot.Width=1.5;
% 
% pard.PhotonRatios.object=struct('Style','edit','String','1 1');%,'Callback',{{@obj.switchvisible,p}});   
% pard.PhotonRatios.position=[6,2.5];
% pard.PhotonRatios.Optional=true;
% pard.PhotonRatios.Width=1.5;


pard.syncParameters={{'cal_3Dfile','',{'String'}}};

pard.plugininfo.type='WorkflowFitter';
pard.plugininfo.description='Global fitter for multiple channels.';
end