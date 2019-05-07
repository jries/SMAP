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
           rrange=500;
           px=2; %nm   
           edge=rrange; %nm
           [locs1,~, hroi]=obj.locData.getloc({'xnm','ynm'},'layer',p.ch1.Value,'Position','roi');
           rxrange=[min(locs1.xnm) max(locs1.xnm)];
           ryrange=[min(locs1.ynm) max(locs1.ynm)];
           
           if p.ch2.Value>1
               locs2=obj.locData.getloc({'xnm','ynm'},'layer',p.ch2.Value-1,'Position','roi');
               rxrange2=[min(locs2.xnm) max(locs2.xnm)];
               ryrange2=[min(locs2.ynm) max(locs2.ynm)];
               rxrange=max(rxrange,rxrange2);
               ryrange=max(ryrange,ryrange2);
           else
               locs2=[];
           end
           
           rx=rxrange(1)-edge:px:rxrange(end)+edge;
           ry=ryrange(1)-edge:px:ryrange(end)+edge;
           
           img1=histcounts2(locs1.xnm,locs1.ynm,rx,ry);
           gr1=fftshift(ifft2(  abs(fft2(img1)).^2));
           
           %mask
           pm=round(([mean(rx) mean(ry)]-p.sr_pos)/px);
           mask=hroi.createMask;
           maskrs=imresize(mask',p.sr_pixrec/px);
           mp=round(size(maskrs)/2);
           rmx=-round(size(img1,1)/2):-round(size(img1,1)/2)+size(img1,1)-1;
           rmy=-round(size(img1,2)/2):-round(size(img1,2)/2)+size(img1,2)-1;
           mcut=maskrs(mp(1)+pm(1)+rmx,mp(2)+pm(2)+rmy);
           mcut(mcut>0)=1;
           mgr=fftshift(ifft2(  abs(fft2(mcut)).^2));
%            density
           A=sum(mcut(:));
           rho=length(locs1.xnm)/A;
           
           %normalize
           gr1n=gr1/rho^2./mgr;
           s=size(gr1n);mp=round(s/2);
           nmax=rrange/px;
           gr1nsm=gr1n(mp(1)-nmax:mp(1)+nmax,mp(2)-nmax:mp(2)+nmax);

           [gr1r,norm]=radialsum(gr1nsm);
           gr1r(end)=[];norm(end)=[];
           nrx=px/2:px:nmax*px+px;
           xaxnm=nrx;
           axpc=obj.initaxis('Pair correlation');
           
           sstart=3;
           gr1rnorm=gr1r(sstart:end)./norm(sstart:end);
           xaxnms=xaxnm(sstart:end);
           if ~isempty(locs2)
                img2=histcounts2(locs2.xnm,locs2.ynm,rx,ry);
                gr2=fftshift(ifft2(  abs(fft2(img2)).^2));
                gr12=fftshift(ifft2(  abs(fft2(img1)).*abs(fft2(img2))));
 
                rho2= length(locs2.xnm)/A;
                gr2n=gr2/rho2^2./mgr;
                gr12n=gr12/rho/rho2./mgr;
                
                gr2nsm=gr2n(mp(1)-nmax:mp(1)+nmax,mp(2)-nmax:mp(2)+nmax);

               [gr2r]=radialsum(gr2nsm);
                gr2r(end)=[];
                gr2rnorm=gr2r(sstart:end)./norm(sstart:end);
                gr12nsm=gr12n(mp(1)-nmax:mp(1)+nmax,mp(2)-nmax:mp(2)+nmax);
                    
               [gr12r]=radialsum(gr12nsm);
                gr12r(end)=[];
                gr12rnorm=gr12r(sstart:end)./norm(sstart:end);
               plot(axpc,xaxnms,gr1rnorm,xaxnms,gr2rnorm,xaxnms,gr12rnorm,xaxnms,1+0*xaxnms,'k')
                legend(axpc,'g1','g2','g12')
           else
                      
            plot(axpc,xaxnms,gr1rnorm,xaxnms,1+0*xaxnms,'k')%,xaxnm,dRn);
           end
           xlabel(axpc,'r (nm)');
           ylabel(axpc,'g(r)');
           
%            Rf=myRipleysK(locs1.xnm,locs1.ynm,nrx,[rx(1) rx(end)  ry(1)  ry(end)]);
           R=myRipleysK(locs1.xnm,locs1.ynm,nrx,[rx(1)+edge rx(end)-edge  ry(1)+edge  ry(end)-edge]);
           L=sqrt(R/pi);
           
           Rpc=cumsum(gr1r)*px^2;
           Lpc=sqrt(Rpc/pi);
           
%            if ~isempty(locs2)
%                
%                R2=myRipleysK(locs2.xnm,locs2.ynm,nrx,[rx(1)+edge rx(end)-edge  ry(1)+edge  ry(end)-edge]);
%                L2=sqrt(R/pi);
%                Rpc2=cumsum(gr2r)*px^2;
%                Lpc2=sqrt(Rpc/pi);
%                       R=myRipleysK(locs1.xnm,locs1.ynm,nrx,[rx(1)+edge rx(end)-edge  ry(1)+edge  ry(end)-edge]);
%            L=sqrt(R/pi);
%            
%            Rpc=cumsum(gr1r)*px^2;
%            Lpc=sqrt(Rpc/pi);
%            else
%            end
          
           ax=obj.initaxis('Ripley');
           plot(ax,nrx,L-nrx',nrx,Lpc-nrx')%,nrx,Lf-Lm)%,nrx,Lpc-nrx')
           legend(ax,'Ripley','Riplex from PC with edge correction')
           ylabel(ax,'L(r)-r');xlabel(ax,'r (nm)');
           
%            ax=obj.initaxis('Ripley from PC');
%            plot(ax,nrx,Lpc-nrx')%,nrx,Lf-Lm)%,nrx,Lpc-nrx')
% %            legend(ax,'Ripley','Riplex from PC with edge correction')
%            ylabel(ax,'L(r)-r');xlabel(ax,'r (nm)');
           
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        

    end
end




function pard=guidef
pard.ch1t.object=struct('String','Ch 1','Style','text');
pard.ch1t.position=[2,1];
pard.ch1t.Width=0.3;
pard.ch1.object=struct('String',{{'layer 1','layer 2','layer 3'}},'Style','popupmenu');
pard.ch1.position=[2,1.3];
pard.ch1.Width=1;
pard.ch2t.object=struct('String','Ch 2','Style','text');
pard.ch2t.position=[2,3];
pard.ch2t.Width=0.3;
pard.ch2.object=struct('String',{{'none','layer 1','layer 2','layer 3'}},'Style','popupmenu');
pard.ch2.position=[2,3.3];
pard.ch2.Width=1;

pard.rranget.object=struct('String','Range r (nm)','Style','text');
pard.rranget.position=[3,1];
pard.rranget.Width=1;
pard.rrange.object=struct('String','1000','Style','edit');
pard.rrange.position=[3,2];
pard.rrange.Width=.5;

pard.rpixt.object=struct('String','Pixelsize (nm)','Style','text');
pard.rpixt.position=[3,3];
pard.rpixt.Width=1;
pard.rpix.object=struct('String','2','Style','edit');
pard.rpix.position=[3,4];
pard.rpix.Width=.5;

pard.plugininfo.name='Spatial Point Patterns';
pard.plugininfo.description= 'Calculates spatial statistics based on pair correlation and Ripleys K function.';
pard.plugininfo.type='ProcessorPlugin';

end