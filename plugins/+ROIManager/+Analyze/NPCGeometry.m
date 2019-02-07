classdef NPCGeometry<interfaces.DialogProcessor&interfaces.SEProcessor
    properties
    end
    methods
        function obj=NPCGeometry(varargin)        
            obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={};
            obj.showresults=true;
        end
        
        function out=run(obj,p)  
            out=runintern(obj,p);
            
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end

    end
end


function out=runintern(obj,p)
out=[];
sites=obj.SE.sites;

% dt=getFieldAsVectorInd(sites,'evaluation.NPCgeomtryQuantify.templatefit.fitted',4);
% ind = dt>20&dt<100;
% sites = sites(ind);
ff='%2.1f';

if isfield(sites(1).evaluation.NPCgeomtryQuantify,'profile') %z-data is there

z0=getFieldAsVector(sites,'evaluation.NPCgeomtryQuantify.profile.Gaussfit.b');
sigma=getFieldAsVector(sites,'evaluation.NPCgeomtryQuantify.profile.Gaussfit.c');
d=getFieldAsVector(sites,'evaluation.NPCgeomtryQuantify.profile.Gaussfit.d');
draw=getFieldAsVector(sites,'evaluation.NPCgeomtryQuantify.profile.Gaussfitraw.d');
dt=getFieldAsVectorInd(sites,'evaluation.NPCgeomtryQuantify.templatefit.fitted',4);

ax0=obj.initaxis('sigma_z');
n=0:1:25;
histogram(sigma,n); xlabel('sigma (nm)')
title(['sigma z: ' num2str(mean(sigma),ff) '\pm' num2str(std(sigma),ff)])
dindin=d>20&d<80;
drawindin=d>20&d<80;
dtindin=d>20&d<80;

d=d(dindin);
draw=draw(drawindin);
dt=dt(dtindin);
n=20:1:80;

ax1=obj.initaxis('d_Gauss');
histogram(abs(d),n); xlabel('d (nm)')
hn=histcounts(abs(d),n);
nf=n(1:end-1)+(n(2)-n(1))/2;
fitp=fit(nf',hn','gauss1');
hold on
plot(nf,fitp(nf))
hold off
title(['d2G: med ' num2str(median(d),ff) ', mean ' num2str(mean(d),ff) '\pm' num2str(std(d),ff) ', fit ' num2str(fitp.b1,ff) '\pm' num2str(fitp.c1/sqrt(2),ff)] )

d2G=fitp.b1;
d2Gerr=fitp.c1/sqrt(2);

ax2=obj.initaxis('d_raw');
histogram(abs(draw),n); xlabel('d (nm)')

hn=histcounts(abs(draw),n);
nf=n(1:end-1)+(n(2)-n(1))/2;
fitp=fit(nf',hn','gauss1');
hold on
plot(nf,fitp(nf))
hold off
title(['draw2G: med ' num2str(median(draw),ff) ', mean ' num2str(mean(draw),ff) '\pm' num2str(std(draw),ff) ', fit ' num2str(fitp.b1,ff) '\pm' num2str(fitp.c1/sqrt(2),ff)] )

draw=fitp.b1;drawerr=fitp.c1/sqrt(2);

ax3=obj.initaxis('d_t');
histogram(abs(dt),n); xlabel('d (nm)')
% title(['template d: median' num2str(median(dt),ff) '\pm' num2str(std(dt),ff)])

hn=histcounts(abs(dt),n);
nf=n(1:end-1)+(n(2)-n(1))/2;
fitp=fit(nf',hn','gauss1');
hold on
plot(nf,fitp(nf))
hold off
title(['dtemplate: med ' num2str(median(dt),ff) ', mean ' num2str(mean(dt),ff) '\pm' num2str(std(dt),ff) ', fit ' num2str(fitp.b1,ff) '\pm' num2str(fitp.c1/sqrt(2),ff)] )

% ind = dt<80&dt>20;
zt=getFieldAsVectorInd(sites,'evaluation.NPCgeomtryQuantify.templatefit.fitted',3);
dpl=d';
indzin=dindin;
ax4=obj.initaxis('d vs z');
hold off
plot(zt(indzin),abs(dpl),'.')
fline=fit(zt(indzin),abs(dpl),'poly1');
hold on
plot(zt(indzin),fline(zt(indzin)),'r')
xlabel('z');ylabel('distance')
title(['d(z=0) fit: ' num2str(fline.p2,ff) ' Corr:' num2str(corr(abs(d'),zt(indzin')))]);
ylim([25 100])
xlim([-200 200])
end

pearsonc=corr(abs(d'),zt(indzin'));

R0=getFieldAsVector(sites,'evaluation.NPCgeomtryQuantify.Rfit');
ax5=obj.initaxis('R');
ff='%2.3f';
% figure(112)
hold off
rn=40:0.5:65;
histogram(abs(R0),rn); xlabel('R (nm)')
title(['fitted radius: ' num2str(mean(R0),ff) '\pm' num2str(std(R0),ff) '\pm' num2str(std(R0)/length(R0),2)])
xlabel('radius (nm)')
xlim([45 60])
ylabel('counts')
hh=histcounts(R0,rn);
fitp=fit(rn(1:end-1)'+(rn(2)-rn(1))/2,hh','gauss1');
hold on
plot(rn,fitp(rn))
% fitp

aca=0;aca1=0;aca2=0;acc12=0;
for k=1:length(sites)
    if isfield(sites(k).evaluation.NPCgeomtryQuantify.angular,'actheta')
    ach=sites(k).evaluation.NPCgeomtryQuantify.angular.actheta;
    aca=ach+aca;
    end
    if isfield(sites(k).evaluation.NPCgeomtryQuantify.angular,'ac1')
        aca1=aca1+sites(k).evaluation.NPCgeomtryQuantify.angular.ac1;
        aca2=aca2+sites(k).evaluation.NPCgeomtryQuantify.angular.ac2;
        acc12=acc12+sites(k).evaluation.NPCgeomtryQuantify.angular.cc12;
    end
end
aca=aca/length(sites);
tn=sites(1).evaluation.NPCgeomtryQuantify.angular.thetan;
axc=obj.initaxis('corr');
aca=aca/2;
% norm=length(tn)-(1:length(tn));
hold off
midp=(length(tn)+1)/2;
dtn=tn(2)-tn(1);
tnn=[tn(midp:end)-2*pi-dtn tn(1:midp-1)];
tnn=tnn/pi*180;
acan=[aca(midp:end) aca(1:midp-1)];
% plot(tn(2:end)*50,aca(2:end)/5);
plot(tnn,acan);
ampm=max(aca(2:end));
ampmin=min(aca);
txtshift=[];
if numel(aca1)>1
    aca1=aca1/length(sites);
    aca2=aca2/length(sites);
    acc12=acc12/length(sites);
    aca2n=[aca2(midp:end) aca2(1:midp-1)];
    aca1n=[aca1(midp:end) aca1(1:midp-1)];
    acc12n=[acc12(midp:end) acc12(1:midp-1)];
    hold on
    plot(tnn,aca1n);
    plot(tnn,aca2n);
    plot(tnn,acc12n);
  
    ampm=max(ampm,max(aca1(2:end)));
    ampm=max(ampm,max(aca2(2:end)));
    ampm=max(ampm,max(acc12(2:end)));
     [ddd,cid]=getcorrangleglobal(tnn,acc12n);
    legend('all','ring1','ring2','cross-corr','cos fit')    
    ff='%2.2f';
     txtshift=[', shift: ' num2str(ddd,ff) '° ± ' num2str((cid(2)-cid(1))/2,ff) '° (95% confidence)'];
      plot([0 0],[ampmin ampm*1.1],'k')
end
ylim([ampmin ampm*1.1]);
xlim([tnn(1) tnn(end)])
ax=gca;
% ax.XTick=-180:45:180;
ax.XAxis.MinorTickValues=-180:45/4:180;
ax.XAxis.TickValues=-180:45:180;
ax.XMinorTick='on';
ax.XGrid='on';
xlabel('angle (theta) ')
ylabel('auto/cross correlation')
title(['angular correlation. Nlocs=' num2str(length(d)) txtshift])
if p.copytopage
    sm=3;sn=3;
    f=figure;
    f.Renderer='painters';
    axt=ax0.copy;
    axt.Parent=f;
    subplot(sm,sn,2,axt)
    
    axt=ax1.copy;
    axt.Parent=f;
    subplot(sm,sn,5,axt) 
    
    axt=axc.copy;
    axt.Parent=f;
   subplot(sm,sn,[7,8,9],axt)
    
    axt=ax3.copy;
    axt.Parent=f;
    subplot(sm,sn,6,axt)
    
    axt=ax4.copy;
    axt.Parent=f;
    subplot(sm,sn,3,axt)
    
    axt=ax5.copy;
    axt.Parent=f;
    subplot(sm,sn,1,axt)
    
    axt=ax2.copy;
    axt.Parent=f;
    subplot(sm,sn,4,axt)
end

%copy to clipboard

filen=obj.SE.files(sites(1).info.filenumber).name;
 lastfile=obj.getPar('lastSMLFile');
 disp([lastfile sprintf('\n') 'N R dR sigma dsigma drawfit drawstdfit d2Gfit d2Gstdfit PearsonC angle dangle95'])
clipboard('copy',[filen sprintf('\t') lastfile sprintf([ '\t' num2str(length(d)) '\t' num2str(mean(R0)) '\t' num2str(std(R0)) ...
    '\t' num2str(mean(sigma)) '\t' num2str(std(sigma)) '\t' num2str(draw) '\t' num2str(drawerr) ...
    '\t' num2str(d2G) '\t' num2str(d2Gerr) '\t' num2str(pearsonc) '\t' num2str(ddd) '\t' num2str((cid(2)-cid(1))/2)])])

end

function [d,cid,fitp]=getcorrangleglobal(ti,cci)
  t=ti;cc=cci;
  mp=round(length(ti)/2);
  win=12;
  t(mp-win:mp+win)=[];
  cc(mp-win:mp+win)=[];
  ft=fittype('(b1+b2*x+b3*x.^2+b4*x.^3)+(a1+a2*x+a3*x.^2+a4*x.^3).*cos(2*pi*(x-d)/45)');
  [~,imcc]=max(cci);
  tm=ti(imcc);
  sp=[(max(cc)-min(cc))/2,0,0,0, min(cc)+1,0,0,0,tm];
  lb=-inf*ones(size(sp));
  lb(1)=0;
  fitp=fit(t',cc',ft,'StartPoint',sp,'Lower',lb);
  t(round(length(t)/2))=nan;
  plot(t,fitp(t))
  d=fitp.d;
  ci=confint(fitp);
  cid=ci(:,end);
%   plot(t,(fitp.b1+fitp.b2*t+fitp.b3*t.^2+fitp.b4*t.^3))
%   plot(t,(fitp.a1+fitp.a2*t+fitp.a3*t.^2+fitp.a4*t.^3+fitp.b1+fitp.b2*t+fitp.b3*t.^2+fitp.b4*t.^3))
%  fitp
end


function pard=guidef(obj)
pard.t1.object=struct('String','Plot results from evaluator: NPCgeometryQuantify','Style','text');
pard.t1.position=[1,1];
pard.t1.Width=4;

pard.copytopage.object=struct('String','Copy to own page','Style','checkbox');
pard.copytopage.position=[2,1];
pard.copytopage.Width=2;


pard.plugininfo.type='ROI_Analyze';

end