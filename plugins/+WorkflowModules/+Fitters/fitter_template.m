classdef fitter_template<interfaces.WorkflowFitter
%     template to integrate a custom fitter in SMAP
    properties
        fitpar %store here parameters for fitting
    end
    methods
        function obj=fitter_template(varargin)
            obj@interfaces.WorkflowFitter(varargin{:})
            obj.inputChannels=1; 
             obj.setInputChannels(1,'frame','segmented ROI stacks');
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function fitinit(obj)
            % this function is called at the beginning once before fitting.
            % Use this to initialiye the fitting parameters and do other
            % time-consuming steps.
            obj.fitpar.roisize=obj.getPar('loc_ROIsize'); % here you can initialize fitting parameters that are used later
            %do other stuff for initialization
            p=obj.getAllParameters; %p also includes all GUI elements
           obj.fitpar.fitpar1=p.fitpar1;
           obj.fitpar.iterations=p.iterations;

            if ~isempty(p.loc_cameraSettings) && p.loc_cameraSettings.EMon %test for EMgain
                obj.fitpar.EMexcessNoise=2;
            else
                obj.fitpar.EMexcessNoise=1;
            end
%             obj.setPar('loc_iterations',p.iterations); % uncomment if one of your parameters is iterations
        end

        function locs=fit(obj,imstack,stackinfo) 
            % this is the function that is called for fitting. it has as
            % input the imstack: stack of single-molecule images, and
            % stackinfo: information on the images, e.g., frame number,
            % position on camera etc.
           
            
            out=fitwrapper(imstack,obj.fitpar); %this function does your fitting
            locs=fit2locs(out,stackinfo,obj.fitpar,imstack); % this function converts the fit output into the SMAP localization data format
            
            if ~isempty(locs)  % copy some information from the stackinfo to the localizations
                fn=fieldnames(stackinfo);
                infonames=setdiff(setdiff(fn,fieldnames(locs)), {'xpix','ypix','frame','Y','X'});
                locs=copyfields(locs,stackinfo,infonames); 
            end
            % the locs structure is then passed on to the next plugin in
            % the workflow
        end

        function initGui(obj) %here you can add functionality to initialize the GUI
            initGui@interfaces.WorkflowFitter(obj);
        end
            
    end
end


function locs=fit2locs(results,stackinfo,fitpar,image)

% here you parse the output of your fitter and convert it to the SMAP
% localizations. You can save as many parameters as you like. But assign at
% least: xpix (x position in camera pixels), ypix, xerrpix (this is the localization precision), yerrpix,
% frame, phot.

if isempty(results)
    locs=[];
    return
end
numl=size(results.P,1);

v1=ones(numl,1,'single');
s=size(image);          
dn=ceil((s(1)-1)/2)*v1;

posx=stackinfo.xpix;
posy=stackinfo.ypix;
frame=stackinfo.frame;
P=results.P;
EMexcess=fitpar.EMexcessNoise;
CRLB=results.CRLB;
LogL=results.LogL;
CRLB(isnan(CRLB))= 0; 
LogL(isnan(LogL))= 0;
CRLB((CRLB)<0)= 0; 
           
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

locs.PSFxpix=P(:,5);
locs.PSFxerr=sqrt(CRLB(:,5));
locs.PSFypix=locs.PSFxpix;

% locs.locpthompson=sqrt((locs.PSFxpix.*locs.PSFypix+1/12*v1)./( locs.phot/EMexcess)+8*pi*(locs.PSFxpix.*locs.PSFypix).^2.* locs.bg./( locs.phot/EMexcess).^2);
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

EMexcess=fitpar.EMexcessNoise; % this is how we deal with the EM excess noise: pretend we have half the photons to mimic the twofold increase in variance
if isempty(EMexcess)
    EMexcess=1;
end
arguments{1}=single(imstack/EMexcess); %image stack
arguments{2}=fitpar.fitpar1; %you can assemble the arguments list of your fitter
arguments{3}=fitpar.iterations;

if 0 %for testing

    [P CRLB LogL]=myfitfunction(arguments{:}); %here you call your fit function. 
else %for testing output zeros
    P=ones(s(3),5,'single');
    P(:,1:2)=s(2)/2+rand(s(3),2); %x,z
    P(:,3)=500; %phot
    CRLB=ones(s(3),5,'single');
    LogL=zeros(s(3),1,'single');

end

out.P=P;
out.CRLB=CRLB;
out.LogL=LogL;
end
        
 

function pard=guidef(obj)

pard.fitpar1t.object=struct('Style','text','String','Par1');
pard.fitpar1t.position=[1,1];
pard.fitpar1.object=struct('Style','edit','String','10');
pard.fitpar1.position=[1,2];
pard.fitpar1.Width=0.5;
pard.fitpar1.TooltipString=sprintf('Fit mode. Fit with constant PSF, free PSF, 3D with astigmatism, asymmetric PSF (for calibrating astigmatic 3D)');

pard.iterationst.object=struct('Style','text','String','iter');
pard.iterationst.position=[1,2.8];
pard.iterations.Width=0.5;
pard.iterations.object=struct('Style','edit','String','50');
pard.iterations.position=[1,3.2];
pard.iterations.TooltipString=sprintf('number of iterations for the GPU fitter (typical 50, use 100-250 for ellipt: PSFx PSFy or 3Dz).');
pard.iterations.Width=0.5;

pard.plugininfo.type='WorkflowFitter';
pard.plugininfo.description='a description';
end