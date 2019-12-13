classdef roi2int_fitG<interfaces.GuiModuleInterface 
    %determines intensity around a localization by a regression of a
    %Gaussian model with fixed positions and sigma. Either amplitude and
    %background or only the amplitude are fitting parameters
    properties
        extension
    end
    methods
        function obj=roi2int_fitG(varargin)
            obj@interfaces.GuiModuleInterface(varargin{:});
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function out=evaluate(obj,p,img,info) %evaluate(obj,roi,bg,dx,dy,PSFxpix,PSFypix)
            out=roi2int_fit_e(p,img,info);
        end
        function prerun(obj,p)
            %TODO include PSF fit
%             global roi2int_fitG_parameters;
%             roi2int_fitG_parameters=obj.getAllParameters;
            warning('off','MATLAB:rankDeficientMatrix');
%             mp=round(sim+1)/2;
%             dn=single(round((roi2int_fitG_parameters.roisize_fit-1)/2));
%             [roi2int_fitG_parameters.X,roi2int_fitG_parameters.Y]=meshgrid(-dn:dn);
            
        end
    end
end


function outp=roi2int_fit_e(p,roi,info)
% global roi2int_fitG_parameters
%weights not implemented? Do htat!
if ~isempty(info.bgim)
    bg=info.bgim;
elseif p.fitonbg
    error('you need to calcualte the background first if you want to subtract it before fitting')
    
%     bg=info.bg;
end
dx=info.dx;
dy=info.dy;
PSFxpix=info.PSFxpix;
PSFypix=info.PSFypix;


sim=size(roi);
if length(sim)==2
    sim(3)=1;
end
mp=round(sim+1)/2;
dn=round((p.roisize_fit-1)/2);
outp=zeros(sim(3),2,'single');
% X=roi2int_fitG_parameters.X;
% Y=roi2int_fitG_parameters.Y;
if p.fixpsf
    PSFxpix=zeros(sim(3),1,'single')+p.psfsize_fit;
    PSFypix=PSFxpix;
end
 nrange=-dn:dn;
if ~p.fitonbg%nargin<7||isempty(bgroi)

    for k=1:sim(3)
        gauss=make2DGaussfast(dx(k),dy(k),PSFxpix(k),PSFypix(k),nrange);
%         gauss=exp((-(dx(k)-X).^2)/2/PSFxpix(k)^2-((dy(k)-Y).^2)/2/PSFypix(k)^2)/pi/PSFxpix(k)/PSFypix(k)/2;
%         weights=sqrt(gauss);
        Xmat=horzcat(gauss(:), gauss(:)*0+1);
        roih=roi(mp(1)-dn:mp(1)+dn,mp(2)-dn:mp(2)+dn,k);
        outp(k,:)=Xmat\roih(:);
        if 0 %p(k,1)>2500 %any(roih(:))%p(k,1)>2500~
            outp(k,:)
            figure(67)
            subplot(2,2,1)
            imagesc(-dn:dn,-dn:dn,roih);
            hold on
            plot(dx(k),dy(k),'+')
            hold off
            subplot(2,2,2);
            imagesc(-dn:dn,-dn:dn,gauss*outp(k,1)+outp(k,2))
            hold on
            plot(dx(k),dy(k),'+')
            hold off
            subplot(2,2,3);
            imagesc(-dn:dn,-dn:dn,gauss*outp(k,1)+outp(k,2)-roih)
            waitforbuttonpress
        end
    end
else %fit bg
%     bgnorm=(2*dn+1)^2;
    nrange=-dn:dn;
    for k=1:sim(3)
        bgh=bg(k);
%         exponent=(-(dx(k)-X).^2)/2/PSFxpix(k)^2-((dy(k)-Y).^2)/2/PSFypix(k)^2;
%         gauss=exp(exponent)/pi/PSFxpix(k)/PSFypix(k)/2;
        gauss=make2DGaussfast(dx(k),dy(k),PSFxpix(k),PSFypix(k),nrange);
%         weights=sqrt(gauss);
        Xmat=horzcat(gauss(:));
%         bgh=bg(mp(1)-dn:mp(1)+dn,mp(2)-dn:mp(2)+dn,k);
        roih=roi(mp(1)-dn:mp(1)+dn,mp(2)-dn:mp(2)+dn,k)-bgh;
        outp(k,1)=Xmat\roih(:);
        outp(k,2)=bgh;
%         p(k,2)=(sum(bgh(:)))/bgnorm;
        if 0%p(k,1)>2500~
            outp(k,:)
            figure(67)
            subplot(2,2,1)
            imagesc(-dn:dn,-dn:dn,roih);
            hold on
            plot(dx(k),dy(k),'+')
            hold off
            subplot(2,2,2);
            imagesc(-dn:dn,-dn:dn,gauss*outp(k,1)+outp(k,2))
            hold on
            plot(dx(k),dy(k),'+')
            hold off
            subplot(2,2,3);
            imagesc(-dn:dn,-dn:dn,gauss*outp(k,1)+outp(k,2)-roih)
            waitforbuttonpress
        end
    end
end
end


function pard=guidef(obj)
pard.t1.object=struct('Style','text','String','roisize');
pard.t1.position=[1,1];
pard.t1.TooltipString='Roi size around localizations for fitting';
% pard.t1.Width=0.5;
pard.roisize_fit.object=struct('Style','edit','String','5');
pard.roisize_fit.position=[1,2];
pard.roisize_fit.TooltipString=pard.t1.TooltipString;
pard.fixpsf.object=struct('Style','checkbox','String','Fix PSF (pix) size to:');
pard.fixpsf.position=[2,1];
pard.fixpsf.Width=3;
pard.fixpsf.TooltipString='fix the PSF during the fit (value in pix). Otherwise the fitted PSF size is used.';

pard.psfsize_fit.object=struct('Style','edit','String','1');
pard.psfsize_fit.position=[2,4];
pard.psfsize_fit.TooltipString=pard.fixpsf.TooltipString;

pard.fitonbg.object=struct('Style','checkbox','String','subtract BG before fitting','Value',0);
pard.fitonbg.position=[3,1];
pard.fitonbg.Width=4;
pard.fitonbg.TooltipString='If selected, the background is subtracted and the fit is performed with an offset=0. Otherwise the background is a fit parameter.';

info.prefix='fit';
info.name='fit';
info.fields={'fit_n','fit_bg'};
pard.plugininfo=info;
pard.plugininfo.type='WorkflowIntensity'; 
pard.plugininfo.description='determines intensity around a localization by a regression of a Gaussian model with fixed positions and sigma. Either amplitude and background or only the amplitude are fitting parameters';
end
