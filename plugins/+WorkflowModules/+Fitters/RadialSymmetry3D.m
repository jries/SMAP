classdef RadialSymmetry3D<interfaces.WorkflowFitter
    properties
        fitpar
    end
    methods
       function obj=RadialSymmetry3D(varargin)
            obj@interfaces.WorkflowFitter(varargin{:})
            obj.inputChannels=2; 
             obj.setInputChannels(2,'frame');
        end
        function pard=guidef(obj)
            pard=guidef;
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
             wx=zeros(s(3),1,'single');
             wy=zeros(s(3),1,'single');
             epsilon=zeros(s(3),1,'single');
             dn=ceil((s(1)-1)/2);
             for k=1:s(3)
                 ims=imstack(:,:,k);
                 bgs=bgstack(:,:,k);
                 [x(k),y(k),epsilon(k),wx(k),wy(k)]=GradientFit3D(ims,dn,dn-1); 
                 bg(k)=mean(bgs(:));
                 phot(k)=sum(ims(:))-bg(k)*s(1)*s(2);
             end
             goodind=abs(x)<3&abs(y)<3;
             
             shiftx=0;%-0.5; %deviation from ground truth
             shifty=0;%-0.5;
             posx=stackinfo.x+shiftx;
             posy=stackinfo.y+shifty;             
            
             locs.frame=stackinfo.frame(goodind);
            locs.xpix=x(goodind)+posx(goodind);
            locs.ypix=-y(goodind)+posy(goodind); %really funny definition 
             locs.PSFxpix=wx(goodind);
             locs.PSFypix=wy(goodind);
             locs.phot=phot(goodind);
             locs.bg=bg(goodind);
%              locs.xerrpix=sqrt((locs.PSFx.*locs.PSFy+1/12)./( locs.phot)+8*pi*(locs.PSFx.*locs.PSFy).^2.* locs.bg./( locs.phot).^2);
              locs.xerrpix=sqrt((locs.PSFxpix.*locs.PSFypix+1/12)./( locs.phot)+8*pi*(locs.PSFxpix.*locs.PSFypix).^2.* locs.bg./( locs.phot).^2);
             locs.peakfindx=posx(goodind);
             locs.peakfindy=posy(goodind);
             locs.gradient3Dellipticity=epsilon(goodind);
             

    end
    end
end


function pard=guidef
pard.plugininfo.type='WorkflowFitter';
pard.plugininfo.description='Adepted from: H. Ma, J. Xu, J. Jin, Y. Gao, L. Lan, and Y. Liu, ?Fast and Precise 3D Fluorophore Localization based on Gradient Fitting.,? Sci. Rep., vol. 5, p. 14335, 2015.';
end