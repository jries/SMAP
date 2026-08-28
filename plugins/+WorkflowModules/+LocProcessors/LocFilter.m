classdef LocFilter<interfaces.WorkflowModule
% Filters localizations before saving according to photons, PSF,
% localization precision, log-likelihood. This can dramatically reduce the
% file size in case a too low cutoff was chosen during the peak finding.
    properties
         pixelsize   
    end
    methods
       function obj=LocFilter(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=1; 
            obj.inputParameters={'loc_ROIsize','loc_iterations'};
%              obj.setInputChannels(1,'frame');
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
        end
        function prerun(obj,p)

           obj.pixelsize=obj.getPar('loc_cameraSettings').cam_pixelsize_um*1000;
            
        end
        function output=run(obj,data,p)
            output=[];
            locs=data.data;%get;
            if data.eof&&isempty(locs)
                output=data;
                return
            elseif isempty(locs)
                disp('no localizations found')
                return
            end
            indin=true(size(locs.frame));
%             if ~isempty(locs)
                %locprec
                if p.check_locprec && isfield(locs,'xerrpix')
                    phot=locs.xerrpix;
                    if isfield(locs,'yerrpix')
                        phot=sqrt((phot.^2+locs.yerrpix.^2)/2);
                    end
                    xerrnm=phot*obj.pixelsize(1);
                    val=p.val_locprec;
                    if length(val)==1
                        val=[0 val];
                    end
                    indinh=xerrnm>=val(1)&xerrnm<=val(2);
                    indin=indin&indinh;
                end
                
                %PSFxnm
                if p.check_psf && isfield(locs,'PSFxpix')
                    psf=locs.PSFxpix;
                    if isfield(locs,'PSFypix')
                        psf=sqrt((psf.^2+locs.PSFypix.^2)/2);
                    end
                    psfnm=psf*obj.pixelsize(1);
                    val=p.val_psf;
                    if length(val)==1
                        val=[obj.pixelsize(1)*.52 val];
                    end
                    indinh=psfnm>=val(1)&psfnm<=val(2);
                    indin=indin&indinh;
                end
                
                %phot
                if p.check_phot && isfield(locs,'phot')
                    phot=locs.phot;
                    val=p.val_phot;
                    if length(val)==1
                        val=[val inf];
                    end
                    indinh=phot>=val(1)&phot<=val(2);
                    indin=indin&indinh;
                end 
                
%                 fn=fieldnames(locs);
%                 for k=1:length(fn)
%                     locsout.(fn{k})=locs.(fn{k})(indin);
%                 end
                %LL
                if p.check_LL && isfield(locs,'logLikelihood')
                    ll=locs.logLikelihood;
                    
                    val=p.val_LL*p.loc_ROIsize^2/2;

                    indinh=ll>=val;
                    indin=indin&indinh;
                end 
                
                if p.check_converged && isfield(locs,'iterations')&&~isempty(p.loc_iterations)
                    indinh=locs.iterations<p.loc_iterations;
                    indin=indin & indinh;
                end
                
                if p.check_convergedxy && isfield(locs,'peakfindx') 
%                  maxfitdist=min(3.5,(p.loc_ROIsize-1)/2);
                    maxfitdist=p.val_convxy;
                    indin=indin & abs(locs.xpix-locs.peakfindx)<maxfitdist & abs(locs.ypix-locs.peakfindy)<maxfitdist;
                end
                
                %remove all locs that contain NaN
                fn=fieldnames(locs);
                for k=1:length(fn)
                    indin(isnan(locs.(fn{k})))=false;
%                     sum(isnan(locs.(fn{k})))
                end    
                
%                 fn=fieldnames(locs);
                for k=1:length(fn)
                    locsout.(fn{k})=locs.(fn{k})(indin);
                end                
                
                output=data;
                output.data=locsout;
%             end
        end

    end
end


function pard=guidef(obj)
pard.txt.object=struct('Style','text','String','Filter ([min max]):');
pard.txt.position=[1,1];
pard.txt.Width=1.3;
pard.txt.Optional=true;
pard.check_locprec.object=struct('Style','checkbox','String','xy-locprec (nm)','Value',0);
pard.check_locprec.position=[2,1];
pard.check_locprec.Width=1.3;
pard.check_locprec.TooltipString=sprintf('Filter localization precision before saving.');
pard.check_locprec.Optional=true;
pard.val_locprec.object=struct('Style','edit','String','100');
pard.val_locprec.position=[2,2.3];
pard.val_locprec.Width=.7;
pard.val_locprec.TooltipString=sprintf('maximum localization precision (nm)');
pard.val_locprec.Optional=true;
pard.check_psf.object=struct('Style','checkbox','String','PSFxy (nm)','Value',0);
pard.check_psf.position=[3,1];
pard.check_psf.Width=1.3;
pard.check_psf.TooltipString=sprintf('Filter size of fitted PSF before saving.');
pard.check_psf.Optional=true;

pard.val_psf.object=struct('Style','edit','String','80 300');
pard.val_psf.position=[3,2.3];
pard.val_psf.Width=.7;
pard.val_psf.TooltipString=sprintf('maximum size of PSF (nm)');
pard.val_psf.Optional=true;

pard.check_phot.object=struct('Style','checkbox','String','Photons','Value',0);
pard.check_phot.position=[4,1];
pard.check_phot.Width=1.3;
pard.check_phot.TooltipString=sprintf('Filter photons before saving.');
pard.check_phot.Optional=true;


pard.val_phot.object=struct('Style','edit','String','200 inf');
pard.val_phot.position=[4,2.3];
pard.val_phot.Width=.7;
pard.val_phot.TooltipString=sprintf('minimum number of photons or vector with minimum and maximum number of photons.');
pard.val_phot.Optional=true;


pard.check_LL.object=struct('Style','checkbox','String','rel. Log Likelihood','Value',0);
pard.check_LL.position=[5,1];
pard.check_LL.Width=1.3;
pard.check_LL.TooltipString=sprintf('Filter log-lieklihood before saving.');
pard.check_LL.Optional=true;

pard.val_LL.object=struct('Style','edit','String','-2');
pard.val_LL.position=[5,2.3];
pard.val_LL.Width=.7;
pard.val_LL.TooltipString=sprintf('Cutoff relative to maximum of log-likelihood distribution (typically 1, not much smaller).');
pard.val_LL.Optional=true;

pard.check_converged.object=struct('Style','checkbox','String','iter< max_iter','Value',0);
pard.check_converged.position=[6,1];
pard.check_converged.Width=1.3;
pard.check_converged.TooltipString=sprintf('Filter fits that did not converge.');
pard.check_converged.Optional=true;

pard.check_convergedxy.object=struct('Style','checkbox','String','|xfit-xpeakfind|<','Value',0);
pard.check_convergedxy.position=[7,1];
pard.check_convergedxy.Width=1.3;
pard.check_convergedxy.TooltipString=sprintf('Filter fits that did not converge to any point close to center of ROI (in pixels)');
pard.check_convergedxy.Optional=true;

pard.val_convxy.object=struct('Style','edit','String','3');
pard.val_convxy.position=[7,2.3];
pard.val_convxy.Width=.7;
pard.val_convxy.TooltipString=sprintf('Filter fits that did not converge to any point close to center of ROI');
pard.val_convxy.Optional=true;



pard.plugininfo.type='WorkflowModule'; 
pard.plugininfo.description='Filters localizations before saving according to photons, PSF, localization precision, log-likelihood. This can dramatically reduce the file size in case a too low cutoff was chosen during the peak finding.';
end