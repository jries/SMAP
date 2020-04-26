classdef zSALM<interfaces.DialogProcessor
%     converts reltive intensities of super- and undercritical channel into
%     z position
    properties
        
    end
    methods
        function obj=zSALM(varargin)   
            obj@interfaces.DialogProcessor(varargin{:}) ;
             obj.showresults=true;
            obj.history=true;
        end
        function out=run(obj,p)
            obj.setPar('undoModule','zSALM');
            notify(obj.P,'backup4undo');
            %get calibration
            fsa=p.assignfield1.selection;
            fua=p.assignfield2.selection;
            fsabg=p.bgfield1.selection;
            fsubg=p.bgfield2.selection;
            locs=obj.locData.getloc({fsa,fua,[fsa 'err'],[fua 'err'],'znm','phot','LLrel','znm_a'},...
                'layer', find(obj.getPar('sr_layerson')),'position','roi');
            if ~isfield(obj.locData.loc,'znm_original')
                obj.locData.setloc('znm_original',obj.locData.loc.znm);
                obj.locData.setloc('locprecznm_original',obj.locData.loc.locprecznm);
            end
%             if ~isfield(obj.locData.loc,'locprecnm_a')
%                  obj.locData.setloc('locprecnm_a',obj.locData.loc.locprecnm);
%             end
            if ~isempty(locs.znm_a)
                locs.znm=locs.znm_a;  
            else
                obj.locData.setloc('znm_a',obj.locData.loc.znm);
                obj.locData.setloc('locprecznm_a',obj.locData.loc.locprecznm);
%                 obj.locData.setloc('locprecnm_a',obj.locData.loc.locprecnm);
            end
            nfit=5e4;
            goodLL=locs.LLrel>-1.3;
            photll=locs.phot(goodLL);
            [~,ind]=sort(photll);
            indref=max(1,length(ind)-nfit);
            photref=photll(ind(indref));
            indbright=(locs.phot>photref) & goodLL & (locs.(fsa)~=0); %take out r=0
            
%             rfactortot=intensitySALM(0)/p.rfactor;
            rfactortot=p.rfactor;
            
            is=locs.(fsa)(indbright);
            iu=locs.(fua)(indbright);
            rsu=double(is./(iu))* rfactortot; %/p.rfactor*intensitySALM(0);
            znm=double(locs.znm(indbright));
            
            dz=10;
            if length(p.zrange)==1
                zrange=-p.zrange:dz:p.zrange;
            elseif length(p.zrange)==2
                zrange=p.zrange(1):dz:p.zrange(2);
            else
                zrange=p.zrange;
            end
                 
            rrange=-0.2:0.01:2.3;
            indf=znm>max(quantile(znm,0.005),zrange(1)) & znm<min(quantile(znm,0.995),zrange(end))...
                & rsu>max(quantile(rsu,0.005),rrange(1)) & rsu<min(quantile(rsu,0.995),rrange(end));
            
            hz=histcounts2(rsu(indf),znm(indf),rrange,zrange);
            
            ax2=obj.initaxis('r vs z');
            hold(ax2,'off')
            h=imagesc(ax2,zrange,rrange,hz);
            axis(ax2, 'xy')
            hold(ax2,'on')
%             ft = fittype('a*exp(-b*x)+c');

%             startp=[1,0.007,0];
%             startp(1)=0.5*exp(startp(2)*median(znm));
%             
%             fitp=fit(znm(indf),rsu(indf),ft,'StartPoint',startp,'Robust','Bisquare','Lower',[0 -inf 0]);
%             
            p.limit=true;
           
            ftsalm=fittype(@(b,x) intensitySALM(x-b,p),'independent','x','coefficients',{'b'});
            startpsalm=quantile(znm(indf),0.1);
            fitpsalm1=fit(znm(indf),rsu(indf),ftsalm,'StartPoint',startpsalm,'Robust','Bisquare');
            %only fit data above coverglass
            zabove = znm>fitpsalm1.b;
            indf=indf&zabove;
            fitpsalm=fit(znm(indf),rsu(indf),ftsalm,'StartPoint',fitpsalm1.b,'Robust','Bisquare');          
            zglass=fitpsalm.b;
            plot(ax2,zrange,fitpsalm(zrange),'r') 
            plot(ax2,[zrange(1), zrange(end)],intensitySALM(0,p)*[1,1],'k') 
            plot(ax2,[zglass zglass],intensitySALM(0,p)*[0, 1],'k')
            title(ax2,['position of glass (nm): ' num2str(zglass,3)]);
            drawnow
           
            isall=obj.locData.loc.(fsa);
            iuall=obj.locData.loc.(fua);
            bgs=obj.locData.loc.(fsabg);
            bgu=obj.locData.loc.(fsubg);
            
            if isfield(obj.locData.loc,[fsa 'err'])
                isallerr=obj.locData.loc.([fsa 'err']);
                iuallerr=obj.locData.loc.([fua 'err']);
            else
                isallerr=sqrt(errN2(isall,bgs));
                iuallerr=sqrt(errN2(iuall,bgu));
            end
            
            
            rall=isall./iuall*rfactortot;%/p.rfactor*intensitySALM(0);
            %checked, this rfactortot should be taken into account
            %correctly during error calculation
            
            %SALM model
            %invert function by spline interpolation    
            %get z coordinates
            dzinterp=0.5;
            zinterp=(-100:dzinterp:2000)';
            rinterp=ftsalm(0,zinterp);
            interpsalm=fit(rinterp,zinterp,'cubicinterp');
            
            zr=interpsalm(rall);
            zmax=zrange(end)+1500;
            zoutofrange=zr>zmax;
            zr(zoutofrange)=zmax; %avoid too large numbers

            obj.locData.loc.znm_SALM=zr;
%             obj.locData.loc.locprecnm=obj.locData.loc.locprecnm_a;
            obj.locData.loc.locprecnm(zoutofrange)=1000; %for grouping.
            obj.locData.loc.znm_a=obj.locData.loc.znm_a-zglass; %correct z astig to put glass to z=0;
            
            %get errors in z
            %calculate dz/dr: use this for r>0.3, here r is a good measure
            %for z
            rdiff=0.01;
            r=(-0.02:rdiff:3)';
            zofr=interpsalm(r);
            dzr=diff(zofr);
            dz_dr=fit(r(1:end-1)+rdiff/2,dzr/rdiff,'cubicinterp');
            rlarge=rall>0.3;
            dzdrv(rlarge,1)=dz_dr(rall(rlarge));
            
            drz=diff(rinterp)/dzinterp;
            dr_dz=fit(zinterp(1:end-1)+dzinterp/2,drz,'cubicinterp');
            dzdrv(~rlarge,1)=1./dr_dz(zr(~rlarge));
            
            %calculate error in z SALM
            dru=-isall./iuall.^2;
            drs=1./iuall;
            dR2=drs.^2.*isallerr.^2+dru.^2.*iuallerr.^2;
            dz2=dzdrv.^2.*dR2;
%             dz2=dz_dr(rall).^2.*dR2;
            indb=isall<1 | iuall<10;
            dz2(zoutofrange|indb)=inf;
            zerrs=sqrt(dz2)*rfactortot;
            errza=obj.locData.loc.locprecznm_a;
            
            znmav=(obj.locData.loc.znm_a./errza.^2+obj.locData.loc.znm_SALM./zerrs.^2)./(1./errza.^2+1./zerrs.^2);
            locprecznmav=1./sqrt(1./errza.^2+1./zerrs.^2);  %divided by two, no idea why, this is not clear
                    
            switch p.fieldznm.Value
                case 1 %weighted average
                    znmnew=znmav;
                    locprecznmnew=locprecznmav;
                case 2 %salm
                    znmnew=obj.locData.loc.znm_SALM;
                    locprecznmnew=zerrs;
                case 3 %astig
                    znmnew=obj.locData.loc.znm_a;
                    locprecznmnew=errza;
            end
            obj.locData.setloc('znm',znmnew);
            obj.locData.setloc('locprecznm',locprecznmnew);
            obj.locData.setloc('locprecznm_SALM',zerrs);
            obj.locData.setloc('znm_asSALM',znmav);
            obj.locData.setloc('locprecznm_asSALM',locprecznmav);
            
            obj.locData.regroup;
            
            %determine maximum position
            rrange=0:0.02:2.3;
            zrange=zrange(1):5:zrange(end);
              hzf=histcounts2(rsu(indf),znm(indf),rrange,zrange);
            
             h=fspecial('gaussian',11,4);
             hf=imfilter(hzf,h);
             factor=5;
             hf=imresize(hf,factor,'cubic');
             [~,linind]=max(hf(:));
             [x,y]=ind2sub(size(hf),linind);
             rmax= x/factor*(rrange(2)-rrange(1))+rrange(1);
             zmax= y/factor*(zrange(2)-zrange(1))+zrange(1);
             plot(ax2,zmax,rmax,'k*')
             disp([rmax,zmax])
             results=sprintf('%3.4f \t',[rmax,zmax]);
             clipboard('copy',results)
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

function fitSALMmodel(z,r,p)
zs=z/1000; %scaling to make r and z similar magnitude

zr=(zs-r)/sqrt(2);
rr=(r+zs)/sqrt(2);

 ftsalm=fittype(@(b,x) (intensitySALM(x*1000-b,p)+x)/sqrt(2),'independent','x','coefficients',{'b'});
            startpsalm=quantile(z,0.1);
            fitpsalm=fit(zr,rr,ftsalm,'StartPoint',startpsalm,'Robust','Bisquare');
end

function loadcall_callback(a,b,obj)
p=obj.getAllParameters;
if isempty(p.cal_3Dfile)
    path=obj.getGlobalSetting('DataDirectory');
    fh=obj.getPar('loc_fileinfo');
    if ~isempty(fh) && ~isempty(fh.imagefile)
        path=fileparts(fh.imagefile);
    end  
    p.cal_3Dfile=[path filesep '*3dcal.mat'];
end
[f,p]=uigetfile(p.cal_3Dfile);
if f
    l=load([p f]);
    if ~isfield(l,'outforfit') && ~isfield(l,'SXY') && ~isfield(l,'cspline')
        msgbox('no 3D data recognized. Select other file.');
    end
    obj.setGuiParameters(struct('cal_3Dfile',[p f]));
    obj.setPar('cal_3Dfile',[p f]);
    
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

pard.rfactort.object=struct('Style','text','String','SA/UA ratio on coverslip');
pard.rfactort.position=[5,1];
pard.rfactort.Width=1.3;
pard.rfactor.object=struct('Style','edit','String','2.16');
pard.rfactor.position=[5,2.3];
pard.rfactor.Width=0.4;

pard.tss.object=struct('String','Field znm','Style','text');
pard.tss.position=[5,3];

pard.fieldznm.object=struct('Style','popupmenu','String',{{'weighted average','SALM','astigmatism'}});
pard.fieldznm.position=[5,3.5];
pard.fieldznm.Width=1.5;
% 
% pard.loadcal.object=struct('Style','pushbutton','String','Load 3D cal','Callback',{{@loadcall_callback,obj}});
% pard.loadcal.position=[6,1];
% pard.loadcal.Width=.75;
% pard.cal_3Dfile.object=struct('Style','edit','String','');
% pard.cal_3Dfile.position=[6,1.75];
% pard.cal_3Dfile.Width=3.25;
% pard.cal_3Dfile.TooltipString=sprintf('3D calibration file for astigmtic 3D. \n Generate from bead stacks with plugin: Analyze/sr3D/CalibrateAstig');


 pard.syncParameters={{'locFields','assignfield1',{'String'}},{'locFields','assignfield2',{'String'}},...
     {'locFields','bgfield1',{'String'}},{'locFields','bgfield2',{'String'}}};%,{'cal_3Dfile','cal_3Dfile',{'String'}}};
            
%             obj.addSynchronization('locFields',[],[],@obj.updateLocFields)
%             obj.addSynchronization('filelist_short',obj.guihandles.dataselect,'String')
pard.plugininfo.type='ProcessorPlugin';

pard.plugininfo.description='converts reltive intensities of super- and undercritical channel into z position.';
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
function err=zerrSALMspline(obj,dz_dr,Ns,Nu,bgs,bgu,zas)
% XXX replace by CRLB from splinePSF
% psf_sa=splinePSF;
% psf_sa.loadmodel(cal_3Dfile,2);
% 
% psf_ua=splinePSF;
% psf_ua.loadmodel(cal_3Dfile,1);

% roi_ua=obj.locData.files.file.savefit.fitparameters.RoiCutterWF.loc_ROIsize;
% errNu2=(psf_ua.crlb(Nu,bgu,zas,roi_ua));

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