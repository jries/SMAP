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
            nfit=5e4;
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
            
            dz=10;
            if length(p.zrange)==1
                zrange=-p.zrange:dz:p.zrange;
            elseif length(p.zrange)==2
                zrange=p.zrange(1):dz:p.zrange(2);
            else
                zrange=p.zrange;
            end
                 
            rrange=-0.2:0.01:2;
            indf=znm>max(quantile(znm,0.005),zrange(1)) & znm<min(quantile(znm,0.995),zrange(end))...
                & rsu>max(quantile(rsu,0.005),rrange(1)) & rsu<min(quantile(rsu,0.995),rrange(end));
    
            
            hz=histcounts2(rsu(indf),znm(indf),rrange,zrange);
            
            ax2=obj.initaxis('r vs z');
            hold(ax2,'off')
            h=imagesc(ax2,zrange,rrange,hz);
            axis(ax2, 'xy')
            hold(ax2,'on')
            ft = fittype('a*exp(-b*x)+c');

            startp=[1,0.007,0];
            startp(1)=0.5*exp(startp(2)*median(znm));
            
            fitp=fit(znm(indf),rsu(indf),ft,'StartPoint',startp,'Robust','Bisquare','Lower',[0 -inf 0]);
            
            
            ftsalm=fittype(@(b,x) intensitySALM(x-b,p),'independent','x','coefficients',{'b'});
            startpsalm=quantile(znm(indf),0.1);
%             ftsalm=fittype(@(b,c,x) intensitySALM(x-b,p)*c,'independent','x','coefficients',{'b','c'});
%             startpsalm=[quantile(znm(indf),0.1) 1];
            
            
            fitpsalm=fit(znm(indf),rsu(indf),ftsalm,'StartPoint',startpsalm,'Robust','Bisquare');
            zglass=fitpsalm.b;
            
            plot(ax2,zrange,fitp(zrange),'m--')  
            plot(ax2,zrange,ftsalm(startpsalm(1),zrange),'y:')
%             plot(ax2,zrange,ftsalm(startpsalm(1),startpsalm(2),zrange),'y:')
            plot(ax2,zrange,fitpsalm(zrange),'r') 
            legend(ax2,'exp','start','salm')
            
            title(ax2,['position of glass (nm): ' num2str(zglass,3)]);
            drawnow

           
            isall=obj.locData.loc.(fsa);
            iuall=obj.locData.loc.(fua);
            rall=isall./iuall;
            
             %exponential model
            zr=-log((rall-fitp.c)/fitp.a)/fitp.b;
           
            
            %SALM model
            %invert function by spline interpolation
            zinterp=(zglass-200:0.5:zrange(end)+500)';
            rinterp=fitpsalm(zinterp);
            interpsalm=fit(rinterp,zinterp,'cubicinterp');
            zr=interpsalm(rall);
            zmax=zrange(end)+500;
            zoutofrange=zr>zmax;
            zr(zoutofrange)=zmax; %avoid too large numbers

            obj.locData.loc.zSALM=zr;
            
            %calculate error of zr
            %use CRLB as error for znm
            % weighted average
            bgs=obj.locData.loc.(fsabg);
            bgu=obj.locData.loc.(fsubg);
%             zerrsexp=zerrSALM(fitp,isall,iuall,bgs,bgu); %now based on exponential fit. Later based on real model?
            
            %calculate dz/dr
            rdiff=0.01;
            r=(-0.02:rdiff:3)';
            zr=interpsalm(r);
            dzr=diff(zr);
            dz_dr=fit(r(1:end-1)+rdiff/2,dzr/rdiff,'cubicinterp');
            
            zerrs=zerrSALMspline(dz_dr,isall,iuall,bgs,bgu);
            zerrs(zoutofrange)=inf;
            errza=obj.locData.loc.locprecznm_a;
            
            %also here change sign of z.
            znmnew=(obj.locData.loc.znm_a./errza+obj.locData.loc.zSALM./zerrs)./(1./errza+1./zerrs);
            locprecznmnew=1./(1./errza+1./zerrs);  %divided by two, no idea why, this is not clear
            obj.locData.setloc('znm',znmnew);
            obj.locData.setloc('locprecznm',locprecznmnew);
            obj.locData.setloc('locprecznm_salm',locprecznmnew);
            obj.locData.regroup;
            out=[];
        end
        
        function initGuiFinal(obj)  
            fields={'assignfield1','assignfield2','bgfield1','bgfield2'};
            fieldnames={'psf_ns','psf_nu','psf_bgs','psf_bgu'};
            for k=1:length(fields)
                stri1=obj.guihandles.(fields{k}).String;
                indg=find(contains(stri1,fieldnames{k}));
                if ~isempty(indg)
                    obj.guihandles.(fields{k}).Value=indg(1);
                end
            end   
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


pard.n1t.object=struct('Style','text','String','n buffer');
pard.n1t.position=[3,1];
pard.n1t.Width=0.7;
pard.n1.object=struct('Style','edit','String','1.33');
pard.n1.position=[3,1.7];
pard.n1.Width=0.4;

pard.n2t.object=struct('Style','text','String','n immersion');
pard.n2t.position=[3,2.2];
pard.n2t.Width=0.7;
pard.n2.object=struct('Style','edit','String','1.78');
pard.n2.position=[3,2.9];
pard.n2.Width=0.4;

pard.lambdat.object=struct('Style','text','String','lambda (nm)');
pard.lambdat.position=[3,3.5];
pard.lambdat.Width=0.7;
pard.lambda.object=struct('Style','edit','String','680');
pard.lambda.position=[3,4.2];
pard.lambda.Width=0.4;

pard.NAt.object=struct('Style','text','String','NA objective');
pard.NAt.position=[4,1];
pard.NAt.Width=0.7;
pard.NA.object=struct('Style','edit','String','1.70');
pard.NA.position=[4,1.7];
pard.NA.Width=0.4;

pard.NAmaskt.object=struct('Style','text','String','NA mask');
pard.NAmaskt.position=[4,2.2];
pard.NAmaskt.Width=0.7;
pard.NAmask.object=struct('Style','edit','String','1.33');
pard.NAmask.position=[4,2.9];
pard.NAmask.Width=0.4;

pard.zranget.object=struct('Style','text','String','Z range (nm)');
pard.zranget.position=[4,3.5];
pard.zranget.Width=0.7;
pard.zrange.object=struct('Style','edit','String','-400 800');
pard.zrange.position=[4,4.2];
pard.zrange.Width=0.8;


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
function err=zerrSALMspline(dz_dr,Ns,Nu,bgs,bgu)
indb=Ns<1 | Nu<10;

% c=max(fitp.c,0);
r=Ns./Nu;
dru=-Ns./Nu.^2;
drs=1./Ns;
dR2=drs.^2.*errN2(Ns,bgs)+dru.^2.*errN2(Nu,bgu);

dz2=dz_dr(r).^2.*dR2;
% dz2=dR2./(-fitp.b.*(r-c)).^2;
% dz2(r<c)=inf;
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