classdef FalconHD2D<interfaces.WorkflowModule
%     FALCON high-density localizer. Speed-optimized version of:	J. Min,
%     C. Vonesch, H. Kirshner, L. Carlini, N. Olivier, S. Holden, S.
%     Manley, J. C. Ye, and M. Unser, ?FALCON: fast and unbiased
%     reconstruction of high-density super-resolution microscopy data.,?
%     Sci. Rep., vol. 4, p. 4577, 2014.';

    properties
        fitpar
        fitfunction
    end
    methods
       function obj=FalconHD2D(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
       end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
       function prerun(obj,p)
           try
            FALCON_GPU_rel3_WF(rand(30,'single'),p.Gsigma,p.speed.selection,false);
            obj.fitfunction=@FALCON_GPU_rel3_WF;
           catch
               obj.fitfunction=@FALCON_CPU_rel3_WF;
               disp('using CPU Falcon fitter');
           end
               
        end
        function dato=run(obj,data,p)
            
            if ~isempty(data.data)
                img=data.data;
                locs=obj.fit(img,p,data.frame);
                dato=data;%.copy;
                dato.data=locs;
            else
                dato=data;
            end
        end
        function locs=fit(obj,imagephot,p,frame)
            debug=false;
      
            [Results,avgImg]= obj.fitfunction(imagephot,p.Gsigma,p.speed.selection,debug);
             shiftx=0;%-0.5; %deviation from ground truth
             shifty=0;%-0.5;
             locs.frame=frame*ones(size(Results,1),1);
             locs.xpix=Results(:,2)+shiftx;
             locs.ypix=Results(:,3)+shifty;
             locs.PSFxpix=p.Gsigma*ones(size(Results,1),1);
             locs.PSFypix=locs.PSFxpix;
             locs.phot=Results(:,4);
             locs.bg=zeros(size(Results,1),1);
             locs.PSF_width_ratio=Results(:,5);
%              locs.xerrpix=sqrt((locs.PSFx.*locs.PSFy+1/12)./( locs.phot)+8*pi*(locs.PSFx.*locs.PSFy).^2.* locs.bg./( locs.phot).^2);
              locs.xerrpix=sqrt((locs.PSFxpix.*locs.PSFypix+1/12)./( locs.phot)+8*pi*(locs.PSFxpix.*locs.PSFypix).^2.* locs.bg./( locs.phot).^2);
             locs.peakfindx=locs.xpix;
             locs.peakfindy=locs.ypix;
    
    end
    end
end

function pard=guidef(obj)
pard.t1.object=struct('Style','text','String','Sigma PSF (cam pix)');
pard.t1.position=[1,1];
pard.Gsigma.object=struct('Style','edit','String','1.2');
pard.Gsigma.position=[1,2];

pard.t2.object=struct('Style','text','String','Speed:');
pard.t2.position=[2,1];
pard.speed.object=struct('Style','popupmenu','String',{{'fast','normal','slow'}});
pard.speed.position=[2,2];
pard.plugininfo.type='WorkflowModule';
pard.plugininfo.description='FALCON high-density localizer. Speed-optimized version of:	J. Min, C. Vonesch, H. Kirshner, L. Carlini, N. Olivier, S. Holden, S. Manley, J. C. Ye, and M. Unser, ?FALCON: fast and unbiased reconstruction of high-density super-resolution microscopy data.,? Sci. Rep., vol. 4, p. 4577, 2014.';
end

