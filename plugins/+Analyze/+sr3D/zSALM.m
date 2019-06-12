classdef zSALM<interfaces.DialogProcessor
    properties
        
    end
    methods
        function obj=zSALM(varargin)   
            obj@interfaces.DialogProcessor(varargin{:}) ;
             obj.showresults=true;

        end
        function out=run(obj,p)
            %get calibration
            fsa=p.assignfield1.selection;
            fua=p.assignfield2.selection;
            fsabg=p.bgfield1.selection;
            fsubg=p.bgfield2.selection;
            locs=obj.locData.getloc({fsa,fua,'znm','phot','LLrel','znm_a'},...
                'layer', find(obj.getPar('sr_layerson')),'position','roi');
            if ~isempty(locs.znm_a)
                locs.znm=locs.znm_a;
            else
                obj.locData.setloc('znm_a',obj.locData.loc.znm);
                obj.locData.setloc('locprecznm_a',obj.locData.loc.locprecznm);
            end
            nfit=3e4;
            goodLL=locs.LLrel>-1.3;
            photll=locs.phot(goodLL);
            [~,ind]=sort(photll);
            indref=max(1,length(ind)-nfit);
            photref=photll(ind(indref));
%             photref=0;
            indbright=(locs.phot>photref)&goodLL;
            
            is=locs.(fsa)(indbright);
            iu=locs.(fua)(indbright);
            rsu=double(is./(iu));
%             rsu=real(log(rsu));
            znm=double(locs.znm(indbright));
            zrange=-0:10:800;
%             rrange=-4.10:0.01:1.5;
            rrange=-0.2:0.01:1.5;
            indf=znm>max(quantile(znm,0.01),zrange(1)) & znm<min(quantile(znm,0.985),zrange(end))...
                & rsu>max(quantile(rsu,0.015),rrange(1)) & rsu<min(quantile(rsu,0.995),rrange(end));
    
            
            hz=histcounts2(rsu(indf),znm(indf),rrange,zrange);
            %make compatible to coordinate system
            hzo=hz';
            ax=obj.initaxis('z vs r');
            hold(ax,'off')
            h=imagesc(ax,rrange,zrange,hzo);
            axis(ax, 'xy')
            hold(ax,'on')
            plot(ax,rrange,0*rrange,'k')
            
                        %
%             try other fit
            
            ax2=obj.initaxis('r vs z');
            hold(ax2,'off')
            h=imagesc(ax2,zrange,rrange,hz);
            axis(ax2, 'xy')
            hold(ax2,'on')
%             ft = fittype('a*exp(-b*x)+c');
            ft = fittype('a*exp(-b*x)+c');
            startp=[.7,.008,0.01];
%             lb=[-inf -inf ];   
%             ft='smoothingspline';
            
            fitp=fit(znm(indf),rsu(indf),ft,'StartPoint',startp,'Robust','Bisquare');
            plot(ax2,zrange,fitp(zrange),'r')  
            plot(ax2,zrange,ft(startp(1),startp(2),startp(3),zrange),'y')
            drawnow
            %later extend to multi-exponential that better describes I(s)
%             ft = fittype('a*exp(-b*(x))');
%             startp=[.45,0.006];
%             lb=[0 0 0];
%                ft = fittype('a*x+b+c*x^2');
%             startp=[-0.005,-1,0];
%             lb=[-inf -inf -inf];   
%              ft = fittype('a*log(x)+b');
%             startp=[-100,-100];
%             lb=[-inf -inf ];   
%             ft = fittype('a*exp(-b*(x))+c*exp(-d*(x))');
%             startp=[.2,0.008, .1 0.001];
%             lb=[0 0 0 0];
%             ft = fittype('a*exp(-b*(x))+c');
%             startp=[.2,0.008, 0];
%             lb=[0 0 -inf];
%              ft = fittype('real(a*log(x)+b+c*x.^2)');
%             startp=[-100,-100,0];
%             lb=[-inf -inf 0];   
%             
%             fitp=fit(rsu(indf),znm(indf),ft,'StartPoint',startp)
%             plot(ax,rrange,fitp(rrange),'r')  
% %             plot(ax,rrange,ft(startp(1),startp(2),rrange),'y')
%             plot(ax,rrange,ft(startp(1),startp(2),startp(3),rrange),'y')
%             drawnow
%               plot(zrange,ft(startp(1),startp(2),startp(3),startp(4),zrange),'y')
            isall=obj.locData.loc.(fsa);
            iuall=obj.locData.loc.(fua);
            rall=isall./iuall;
            
          
            zr=-log((rall-fitp.c)/fitp.a)/fitp.b;
            
%             zr=fitp(rall);
            zr(isnan(zr))=zrange(end)-100;
%             zr=real(log(rall/fitp.a)/(-fitp.b));
%             zr(rsu<0)=1000;
            obj.locData.loc.zSALM=zr;
            
            %calculate error of zr
            %use CRLB as error for znm
            % weighted average
            bgs=obj.locData.loc.(fsabg);
            bgu=obj.locData.loc.(fsubg);
            zerrs=zerrSALM(fitp,isall,iuall,bgs,bgu);
            errza=obj.locData.loc.locprecznm_a;
            %also here change sign of z.
            znmnew=(obj.locData.loc.znm_a./errza+obj.locData.loc.zSALM./zerrs)./(1./errza+1./zerrs);
            locprecznmnew=1./(1./errza+1./zerrs);  %divided by two, no idea why, this is not clear
            obj.locData.setloc('znm',znmnew);
            obj.locData.setloc('locprecznm',locprecznmnew);
            obj.locData.regroup;


%             plot(ax,rrange,0*rrange,'k')
            out=[];
        end
        
        function initGui(obj)           
        end
 
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function pard=guidef(obj)


pard.t1.object=struct('String','Intensity SA','Style','text');
pard.t1.position=[1,1];
pard.t2.object=struct('String','Intensity UA','Style','text');
pard.t2.position=[2,1];

pard.assignfield1.object=struct('Style','popupmenu','String','');
pard.assignfield1.position=[1,2];
pard.assignfield2.object=struct('Style','popupmenu','String','');
pard.assignfield2.position=[2,2];

pard.t3.object=struct('String','bg SA','Style','text');
pard.t3.position=[1,3.5];
pard.t4.object=struct('String','bg UA','Style','text');
pard.t4.position=[2,3.5];

pard.bgfield1.object=struct('Style','popupmenu','String','');
pard.bgfield1.position=[1,4];
pard.bgfield2.object=struct('Style','popupmenu','String','');
pard.bgfield2.position=[2,4];


 pard.syncParameters={{'locFields','assignfield1',{'String'}},{'locFields','assignfield2',{'String'}},...
     {'locFields','bgfield1',{'String'}},{'locFields','bgfield2',{'String'}}};
            
%             obj.addSynchronization('locFields',[],[],@obj.updateLocFields)
%             obj.addSynchronization('filelist_short',obj.guihandles.dataselect,'String')
pard.plugininfo.type='ProcessorPlugin';


end

function err=zerrSALM(fitp,Ns,Nu,bgs,bgu)
indb=Ns<1 | Nu<10;
% Ns(Ns<1)=1;
% Nu(Nu<1)=1;
c=max(fitp.c,0);
r=Ns./Nu;
dru=-Ns./Nu.^2;
drs=1./Ns;
dR2=drs.^2.*errN2(Ns,bgs)+dru.^2.*errN2(Nu,bgu);
dz2=dR2./(-fitp.b.*(r-c)).^2;
dz2(r<c)=inf;
dz2(indb)=inf;
% dz2=(a./r).^2.*dR2;
% dz2=dR2/b^2./r.^2;
err=sqrt(dz2);
end

function dN2=errN2(N,bg)
bg(bg<0)=0;
s_a=1; %sigmapsf/pixelsize
tau=2*pi*(bg)*(s_a^2+1/12)./N;
dN2=N.*(1+4*tau+sqrt(tau./(14*(1+2*tau))));
end