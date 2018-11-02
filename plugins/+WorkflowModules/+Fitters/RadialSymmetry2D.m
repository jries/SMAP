classdef RadialSymmetry2D<interfaces.WorkflowFitter
    properties
        fitpar
    end
    methods
       function obj=RadialSymmetry2D(varargin)
            obj@interfaces.WorkflowFitter(varargin{:})
            obj.inputChannels=2; 
             obj.setInputChannels(2,'frame');
        end
        function pard=guidef(obj)
            pard.plugininfo.type='WorkflowFitter';
            pard.plugininfo.description='Fast CPU-based radial symmetry localizer (2D) according to: R. Parthasarathy, ?Rapid, accurate particle tracking by calculation of radial symmetry centers,? Nat Methods, Jun. 2012.';
          
        end
        function fitinit(obj)
            obj.numberInBlock=0;
        end
        function locs=fit(obj,imstack,bgstack,stackinfo)
            
             s=size(imstack);
             if length(s)==2 
                 s(3)=1;
             end
             if s(3)==0
                 locs=[];
                 return
             end
             
             x=zeros(s(3),1,'single');
             y=zeros(s(3),1,'single');
             bg=zeros(s(3),1,'single');
             phot=zeros(s(3),1,'single');
             sigma=zeros(s(3),1,'single');
             for k=1:s(3)
                 ims=imstack(:,:,k);
                 bgs=bgstack(:,:,k);
                 [x(k),y(k),sigma(k)]=radialcenter(ims); 
                 bg(k)=mean(bgs(:));
                 phot(k)=sum(ims(:))-bg(k)*s(1)*s(2);
             end
             
             dn=ceil((s(1)-1)/2);
             shiftx=0;%-0.5; %deviation from ground truth
             shifty=0;%-0.5;
             posx=stackinfo.x+shiftx;
             posy=stackinfo.y+shifty;
             locs.frame=stackinfo.frame;
             locs.xpix=x-dn+posx-1;
             locs.ypix=y-dn+posy-1;
             locs.PSFxpix=sigma;
             locs.PSFypix=locs.PSFxpix;
             locs.phot=phot;
             locs.bg=bg;
%              locs.xerrpix=sqrt((locs.PSFx.*locs.PSFy+1/12)./( locs.phot)+8*pi*(locs.PSFx.*locs.PSFy).^2.* locs.bg./( locs.phot).^2);
              locs.xerrpix=sqrt((locs.PSFxpix.*locs.PSFypix+1/12)./( locs.phot)+8*pi*(locs.PSFxpix.*locs.PSFypix).^2.* locs.bg./( locs.phot).^2);
             locs.peakfindx=posx;
             locs.peakfindy=posy;
    
    end
    end
end

