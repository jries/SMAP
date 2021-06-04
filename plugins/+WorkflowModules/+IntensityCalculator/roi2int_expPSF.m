classdef roi2int_expPSF<interfaces.GuiModuleInterface 
%    Determines intensity around a localization by a) regression of an
%    experimental PSF model (amplitude and background or only the amplitude
%    free fitting), or b) by multipliation with an experimental PSF model.
    properties
        splinecoeff
        spline
        extension
%         transform
%         p
        mirror
    end
    methods
        function obj=roi2int_expPSF(varargin)
            obj@interfaces.GuiModuleInterface(varargin{:});
        end
        function initGui(obj)
            obj.makeinfobutton;
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function out=evaluate(obj,p,roi,loc)
            out=roi2int_fit_e(obj,p,roi,loc);
        end
        function load3Dfile(obj,a,b)
            calf=obj.getSingleGuiParameter('cal_3Dfile');
            if isempty(calf) %nothing selected
                mainf=obj.getPar('mainfile');
                if ~isempty(mainf)
                    pf=fileparts(mainf);
                    calf=[pf,filesep,'*_3dcal.mat'];
                end
            end
            [file,pfad]=uigetfile(calf);
            if file
                obj.setGuiParameters(struct('cal_3Dfile',[pfad file]));
            end
            load_spline(obj)
        end
        function load_spline(obj)
            f=obj.getSingleGuiParameter('cal_3Dfile');
            if isempty(f)
                return
            end
            obj.spline=load(f);
%             sppos=min(obj.guihandles.splinefields.Value,length(obj.spline.SXY));
            ic=obj.getPar('intensity_channel');
            if isempty(ic)
                ic=obj.extension;
                if isempty(ic)
                    ic='t';
                end
%                 return;
            end
            switch ic
                case {'t' ,'s'}
                    sppos=2;
                case {'r','u'}
                    sppos=1;
            end
%                sppos        
            sind=min(length(obj.spline.SXY(sppos).cspline.coeff),sppos); %did not fix calibrator yet..
            obj.splinecoeff=single(obj.spline.SXY(sppos).cspline.coeff{sind});
            
        end
        function prerun(obj,p)
%              fit3ddir=strrep(pwd,'SMAP','fit3D');
%             if exist(fit3ddir,'file') && ~isdeployed
%             addpath(fit3ddir);
%             end
%             
%             fit3ddir=strrep(pwd,'SMAP','fit3Dcspline');
%             if exist(fit3ddir,'file')&& ~isdeployed
%             addpath(fit3ddir);
%             end
%             obj.p=obj.getGuiParameters;
            obj.load_spline;
            warning('off','MATLAB:lscov:RankDefDesignMat')

%             dn=round((obj.p.roisize_fit-1)/2);
%             zmp=obj.spline.SXY(1).cspline.z0;
%             norm=evalSpline(obj.p.roisize_fit,obj.splinecoeff,1,0,[dn dn zmp]);
%             sppos=obj.guihandles.splinefields.Value;
%             transform=obj.spline.transformation;
%             if isa(transform,'interfaces.LocTransformN')
%                 obj.mirror=transformation.mirror(sppos);
%             else
                
%                 obj.mirror=false(1,2);
%                  mm=transform.tinfo.mirror.targetmirror;
%                  if contains(mm,'up') && sppos==2
%                      obj.mirror(2)=true;
%                  end
%                  if contains(mm,'right') && sppos==2
%                       obj.mirror(1)=true;
%                  end
%             end
%             roi2int_fit_e();
            %TODO include PSF fit
%             global roi2int_fitG_parameters;
%             roi2int_fitG_parameters=obj.getAllParameters;
%             mp=round(sim+1)/2;
%             dn=single(round((roi2int_fitG_parameters.roisize_fit-1)/2));
%             [roi2int_fitG_parameters.X,roi2int_fitG_parameters.Y]=meshgrid(-dn:dn);
            
        end
    end
end


function pout=roi2int_fit_e(obj,p,roi,loc)
global debugon
% persistent splinehere
% if nargin==0
%     splinehere=[];
% end
% if isempty(splinehere)
% end

%weights not implemented? Do htat!
sim=size(roi);
if length(sim)==2
    sim(3)=1;
end

multiply=strcmp(p.fitmode.selection,'multiply');
% if obj.mirror(2)
%     roi(:,:,:)=roi(end:-1:1,:,:);
%     loc.dy=-loc.dy;
% end
% if obj.mirror(1)
%     roi(:,:,:)=roi(:,end:-1:1,:);
%     loc.dx=-loc.dx;
% end
mp=round(sim+1)/2;
dn=round((p.roisize_fit-1)/2);
pout=zeros(sim(3),2,'single');
% n=-dn:dn;
% [X,Y]=meshgrid(n);
% Z=0*X;
% sppos=obj.guihandles.splinefields.Value;
% splinecoeff=obj.spline.SXY(sppos).cspline.coeff{sppos};
zmp=obj.spline.SXY(1).cspline.z0;
dz=obj.spline.SXY(1).cspline.dz;
    if p.fixz0
        z=zeros(sim(3),1)+p.z0;
    else
        z=loc.znm/dz;
    end
    cor=horzcat(loc.dx+dn,loc.dy+dn,-z+zmp);
   
   % template = evalSpline(obj.p.roisize_fit,obj.splinecoeff,1,0,cor);
    template = simSplinePSF_call(p.roisize_fit,obj.splinecoeff,1,0,single(cor));
   numrois=1;
for k=1:sim(3)
    templateh=template(:,:,k);
    roih=roi(mp(1)-dn:mp(1)+dn,mp(2)-dn:mp(2)+dn,k);
    switch p.fitbg.selection
        case 'BG free fit'
            %fit bg
            if multiply
                error('free fit of BG not possible with multiplication mode')
            end
            Xmat=horzcat(templateh(:),ones((2*dn+1)^2,1));          
            if p.normalizeimage
                ph=lscov(Xmat,roih(:),1./sqrt(roih(:)));
                pout(k,:)=ph;
            else
                pout(k,:)=Xmat\roih(:);
            end
            
        otherwise
            if isfield(loc,'numrois') %grouped data
                numrois=loc.numrois(k);
            end
            switch p.fitbg.selection
                case 'BG from localizations'
                    bg=loc.bg(k);
                case 'BG fixed to:'
                    bg=p.bgsetvalue*numrois;
                case 'BG calculated above'
                    if isempty(loc.bgim)
                        error('you need to select calcualte BG above')
                    else
                    bg=loc.bgim(k)*numrois;
                    end
                otherwise
                    warning('not implemented')      
            end
        %bg fixed
        if multiply
            if p.normalizeimage
                weights=1./sqrt(roih(:));
                pout(k,1)=sum((roih(:)-bg).*templateh(:).*weights)/mean(weights)/sum(templateh(:));
            else
                pout(k,1)=sum((roih(:)-bg).*templateh(:))/sum(templateh(:));
            end
        else
        Xmat=templateh(:);
        if p.normalizeimage
            ph=lscov(Xmat,(roih(:)-bg),1./sqrt(roih(:)));
            pout(k,1)=ph;
        else
            pout(k,1)=Xmat\(roih(:)-bg);
        end
        end
        pout(k,2)=bg;
        
    end
end

if ~isempty(debugon) && debugon%loc.frame(1)==60
    
    roia=roi(mp(1)-dn:mp(1)+dn,mp(2)-dn:mp(2)+dn,:);
    tp=template;
    for k=1:size(roi,3)
        tp(:,:,k)=tp(:,:,k)*pout(k,1)+pout(k,2);
    end
    impl=[roia,tp;roia-tp,0*tp];
    f=figure(88);
    imx(impl,'Parent',f)
    waitforbuttonpress
end
end


function pard=guidef(obj)
pard.t1.object=struct('Style','text','String','roisize');
pard.t1.position=[1,1];
pard.t1.TooltipString='Roi size around localizations for fitting';
pard.t1.Width=0.6;
% pard.t1.Width=0.5;
pard.roisize_fit.object=struct('Style','edit','String','9');
pard.roisize_fit.position=[1,1.6];
pard.roisize_fit.TooltipString=pard.t1.TooltipString;
pard.roisize_fit.Width=0.4;
% 
% pard.t2.object=struct('Style','text','String','spline cal. SXY index:');
% pard.t2.position=[2,1];
% pard.t2.Width=2;
% pard.t2.TooltipString='Roi size around localizations for fitting';
% pard.splinefields.object=struct('Style','popupmenu','String',{{1,2}});
% pard.splinefields.position=[2,3];


pard.loadcal_3Dfile.object=struct('Style','pushbutton','String','load','Callback',@obj.load3Dfile);
pard.loadcal_3Dfile.position=[2,4];
pard.loadcal_3Dfile.Width=1;
pard.loadcal_3Dfile.TooltipString=sprintf('3D calibration file for astigmtic 3D. \n Generate from bead stacks with plugin: Analyze/sr3D/CalibrateAstig');


pard.cal_3Dfile.object=struct('Style','edit','String','');
pard.cal_3Dfile.position=[2,1];
pard.cal_3Dfile.Width=3;
pard.cal_3Dfile.TooltipString=sprintf('3D calibration file for astigmtic 3D. \n Generate from bead stacks with plugin: Analyze/sr3D/CalibrateAstig');

pard.normalizeimage.object=struct('Style','checkbox','String','normalize image to Poisson noise');
pard.normalizeimage.position=[3,1];
pard.normalizeimage.Width=3;
% pard.normalizeimage.TooltipString=sprintf('3D calibration file for astigmtic 3D. \n Generate from bead stacks with plugin: Analyze/sr3D/CalibrateAstig');



pard.fixz0.object=struct('Style','checkbox','String','Fix z pos to frame:');
pard.fixz0.position=[4,1];
pard.fixz0.Width=1.6;
% pard.fixz0.TooltipString='fix the PSF during the fit (value in nm). Otherwise the fitted PSF size is used.';

pard.z0.object=struct('Style','edit','String','0');
pard.z0.position=[4,2.6];
pard.z0.Width=0.4;
% pard.z0.TooltipString=pard.fixpsf.TooltipString;

p(1).value=1:3; p(1).on={}; p(1).off={'bgsetvalue'};
p(2).value=4; p(2).on={'bgsetvalue'}; p(2).off={};
            
pard.fitbg.object=struct('Style','popupmenu','String',{{'BG free fit','BG from localizations','BG calculated above','BG fixed to:'}},'Value',1,...
    'Callback',{{@obj.switchvisible,p}});
pard.fitbg.position=[1,2];
pard.fitbg.Width=2;

pard.fitmode.object=struct('Style','popupmenu','String',{{'fit','multiply'}},'Value',1);
pard.fitmode.position=[3,4];
pard.fitmode.Width=1;

pard.bgsetvalue.object=struct('Style','edit','String','0','Visible','off');
pard.bgsetvalue.position=[1,4.5];
pard.bgsetvalue.Width=0.5;

% pard.fitbg.TooltipString='If selected, the background is subtracted and the fit is performed with an offset=0. Otherwise the background is a fit parameter.';

info.prefix='fit';
info.name='fit';
info.fields={'psf_n','psf_bg'};
pard.plugininfo=info;
pard.plugininfo.type='WorkflowIntensity'; 
pard.syncParameters={{'cal_3Dfile','cal_3Dfile',{'String'}}};
pard.plugininfo.description='Determines intensity around a localization by a) regression of an experimental PSF model (amplitude and background or only the amplitude free fitting), or b) by multipliation with an experimental PSF model.';
end
