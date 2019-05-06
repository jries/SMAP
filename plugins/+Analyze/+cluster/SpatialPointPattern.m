classdef SpatialPointPattern<interfaces.DialogProcessor
 
    properties
        
    end
    methods
        function obj=SpatialPointPattern(varargin)    
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.inputParameters={'sr_pixrec','numberOfLayers','sr_pos','sr_size','layers','sr_layerson'};
%             obj.history=true;    
            obj.showresults=true;
        end
        
        function out=run(obj,p)
            %include tilint?
           out=[];
           layers=[1 2];
           [locs1,~, hroi]=obj.locData.getloc({'xnm','ynm'},'layer',layers(1),'Position','roi');
           edge=1000; %nm
           px=5; %nm
           rx=min(locs1.xnm)-edge:px:max(locs1.xnm)+edge;
           ry=min(locs1.ynm)-edge:px:max(locs1.ynm)+edge;
           img1=histcounts2(locs1.xnm,locs1.ynm,rx,ry);
           
           gr1=fftshift(ifft2(  abs(fft2(img1)).^2));
           mask=hroi.createMask;
           maskrs=imresize(mask,p.sr_pixrec/px);
           mgr=fftshift(ifft2(  abs(fft2(maskrs,size(img1,1),size(img1,2))).^2));
           A=sum(maskrs(:));
           rho=length(locs1.xnm)/A;
%            rho=1;
           gr1n=gr1/rho^2./mgr;
           s=size(gr1n);mp=round(s/2);
           nmax=100;
           gr1nsm=gr1n(mp(1)-nmax:mp(1)+nmax,mp(2)-nmax:mp(2)+nmax);
           
            
           [gr1r,norm]=radialsum(gr1nsm);
           xaxnm=linspace(1,nmax*px,length(gr1r));
           axpc=obj.initaxis('Pair correlation');plot(axpc,xaxnm,gr1r./norm);
           
           nrx=0:px:nmax*px;
           R=myRipleysK(locs1.xnm,locs1.ynm,nrx,[rx(1)+edge rx(end)-edge ry(1)+edge ry(end)-edge]);
           
           mgr1r=mgr(mp(1)-nmax:mp(1)+nmax,mp(2)-nmax:mp(2)+nmax);
           [mask1r,norm]=radialsum(mgr1r);
           
           mRn=mask1r(1:length(R))./norm(1:length(R))/A;
           Rn=R.*mRn;
           
           L=sqrt(R/pi);
           ax=obj.initaxis('Ripley');
           plot(ax,nrx,L-nrx')
           ylabel(ax,'L(r)-r');xlabel(ax,'r (nm)');
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        

    end
end




function pard=guidef
pard.countingregion.object=struct('String','Gauss wighted counting |circle counting','Style','popupmenu','Value',2);
pard.countingregion.object.TooltipString=sprintf('count Gauss-weighted locs (more accurate) or locs in circle/cylinder (faster)');
pard.countingregion.position=[1,1];
pard.countingregion.Width=2;

pard.countwhat.object=struct('String','all|layers','Style','popupmenu','Value',2);
pard.countwhat.object.TooltipString=sprintf('count all localizations or per layer (use visible ones as reference)');
pard.countwhat.position=[2,1];
pard.countwhat.Width=1;

pard.allfiles.object=struct('String','on all files','Style','checkbox','Value',0);
pard.allfiles.object.TooltipString=sprintf('Apply on all files together');
pard.allfiles.position=[2,2];
pard.allfiles.Width=1;

pard.texta.object=struct('String','size in x,y (nm)','Style','text');
pard.texta.position=[3,1];
pard.countingsize_xy.object=struct('String','12','Style','edit');
pard.countingsize_xy.position=[3,2];
pard.countingsize_xy.isnumeric=1;
pard.countingsize_xy.object.TooltipString=sprintf('radius of circle or sigma of gauss in lateral direction');


pard.text1.object=struct('String','size in z (nm)','Style','text');
pard.text1.position=[4,1];
pard.countingsize_z.object=struct('String','24','Style','edit');
pard.countingsize_z.position=[4,2];
pard.countingsize_z.isnumeric=1;
pard.countingsize_z.object.TooltipString=sprintf('size of cylinder or sigma of gauss in z direction');
pard.plugininfo.name='Spatial Point Patterns';
pard.plugininfo.description= 'density_calculator looks at the neighborhood and counts number of neighbours. If grouped or ungrouped data is used depends on setting in layers.';
pard.plugininfo.type='ProcessorPlugin';

end