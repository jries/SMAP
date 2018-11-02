classdef NPCLabelingQuantify<interfaces.SEEvaluationProcessor
    properties
        savedevals
    end
    methods
        function obj=NPCLabelingQuantify(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function makeGui(obj)
            makeGui@interfaces.SEEvaluationProcessor(obj);
%             obj.guihandles.saveimagesb.Callback={@saveimagesb_callback,obj};
        end
        function out=run(obj,p)
%             try
            out=runintern(obj,p);
%             catch err
%                 err
%                     out.numcornersfiltered=0;
%                 out.numcorners=0;
%                 out.numfoundint=0;
%                 out.numfoundrat=0;
%                 out.numbercornerassined=0;
%                 out.rotation=0;
%             end
         
        end
        
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end

function pard=guidef
pard.Rt.object=struct('Style','text','String','R (nm):');
pard.Rt.position=[1,1];
pard.Rt.Width=1;

pard.R.object=struct('Style','edit','String','50');
pard.R.position=[1,2];
pard.R.Width=1;

pard.dRt.object=struct('Style','text','String','dR:');
pard.dRt.position=[1,3];
pard.dRt.Width=1;

pard.dR.object=struct('Style','edit','String','20');
pard.dR.position=[1,4];
pard.dR.Width=1;
% pard.dualchannel.object=struct('Style','checkbox','String','dual channel','Value',0);
% pard.dualchannel.position=[2,1];
% pard.dualchannel.Width=2;
% 
% 
% pard.sqrtfit.object=struct('Style','checkbox','String','fit sqrt(img)','Value',1);
% pard.sqrtfit.position=[3,1];
% pard.sqrtfit.Width=2;
% 
% 
% pard.t2.object=struct('Style','text','String','gaussfac for imfit (l1, l2)');
% pard.t2.position=[4,1];
% pard.t2.Width=3;
% 
% pard.gaussfac_imfit.object=struct('Style','edit','String','1, 1');
% pard.gaussfac_imfit.position=[4,4];
% pard.gaussfac_imfit.Width=1;
% 
% pard.fit_sigma.object=struct('Style','checkbox','String','Fit sigma of ring','Value',0);
% pard.fit_sigma.position=[5,1];
% pard.fit_sigma.Width=2;

pard.plugininfo.type='ROI_Evaluate';
pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi','layer1_','layer2_','se_sitepixelsize'};
end


function out=runintern(obj,p)
R=p.R;
dR=p.dR;

locs=obj.getLocs({'xnm','ynm','xnm_gt','ynm_gt','locprecnm','frame'},'layer',1,'size',p.se_siteroi(1)/2);
if isempty(locs.xnm)
    out=[];
    return
end


[x0,y0]=fitposring(locs.xnm,locs.ynm,R);
xm=locs.xnm-x0;
ym=locs.ynm-y0;

[xfr,yfr,Rfr]=fitposring(locs.xnm,locs.ynm,[],[x0,y0,R]);
out.radius=Rfr;

step=2*pi/8;
[tha,rhoa]=cart2pol(xm,ym);

minlp=step*R*.4;
inr=rhoa>R-dR&rhoa<R+dR;

inr=inr&locs.locprecnm<minlp;
% th=tha(inr);rho=rhoa(inr);
out.numlocsf=sum(inr);
out.numlocsR=sum(rhoa<R+dR);
out.coordinates=struct('rho',rhoa,'theta',tha,'drho',locs.locprecnm, 'dtheta',locs.locprecnm./rhoa,'x',locs.xnm,'y',locs.ynm,'frame',locs.frame);


if obj.display
   axdat=obj.setoutput('Data');
   tabdat=uitabgroup(axdat.Parent);
   ax1=axes(uitab(tabdat,'Title','assigned'));
   ax1b=axes(uitab(tabdat,'Title','a_t'));
   ax2=axes(uitab(tabdat,'Title','spacing'));
%    ax1= obj.setoutput('assigned');
%    ax1b= obj.setoutput('a_t');
%    ax2= obj.setoutput('spacing');
   if ~isempty(locs.xnm_gt)
        axgt=obj.setoutput('Ground_truth');
       tabdat=uitabgroup(axgt.Parent);
       ax3=axes(uitab(tabdat,'Title','assigned'));
       ax3b=axes(uitab(tabdat,'Title','a_t'));
       ax4=axes(uitab(tabdat,'Title','spacing'));
      
   end
else
    ax1=[];ax2=[]; ax3=[];ax4=[]; ax3b=[]; ax1b=[];
end
  
locsall=obj.locData.getloc({'frame','filenumber'},'layer',1,'removefilter','filenumber','position','all');

frames=locsall.frame(locsall.filenumber==obj.site.info.filenumber);
frames0=frames(frames>0);
numpoints=10;
qq=linspace(0,1,numpoints+1);
timepoints=myquantile(frames0,qq);  % now edges are defined: think about time windows
% maxf=max();
% minf=min(locsall.frame(locsall.filenumber==obj.site.info.filenumber));
% 
% 
% 
% qq=timepoints(1,:)/timepoints(1,end);

    

% timepoints=linspace(minf,maxf,numpoints+1);
[numbercornerassigned,mdt,assigneddirect,timing]=assigntocornersdirect(out.coordinates,inr,8,ax1,timepoints,ax1b);
[numfound, numfound2]=countspacing(out.coordinates,inr,8,ax2);


%when considering maximum gap, also look at next corners
%correlation: when uncertain and moves in one direction, distance to othr
%gets larger
% figure(188);histogram(th,-pi:pi/8:pi)

% histogram(dth/step,30)
% subplot(3,1,2)
% plot(dth,dstep,'+')



out.numfoundint=numfound;
out.numfoundrat=numfound2;
out.numbercornerassigned=numbercornerassigned;
out.numbercornerassigneddirect=assigneddirect;
out.timing=timing;
out.rotation=(mdt);

if ~isempty(locs.xnm_gt)
    [x0gt,y0gt]=fitposring(locs.xnm_gt,locs.ynm_gt,R);
    xmgt=locs.xnm_gt-x0gt;
    ymgt=locs.ynm_gt-y0gt;
    [thagt,rhoagt]=cart2pol(xmgt,ymgt);
    coordinates_gt=struct('rho',rhoagt,'theta',thagt,'drho',locs.locprecnm*0+5, 'dtheta',locs.locprecnm*0+pi/64,'x',locs.xnm_gt,'y',locs.ynm_gt,'frame',locs.frame);
    [numbercornerassigned_gt,~,direct_gt]=assigntocornersdirect(coordinates_gt,inr,8,ax3,timepoints);
    [numgap_gt, numgapf_gt]=countspacing(coordinates_gt,inr,8,ax4);

    locsall=obj.getLocs({'xnm','ynm','xnm_gt','ynm_gt','locprecnm','frame'},'size',p.se_siteroi);
    xmgta=locsall.xnm_gt-x0gt;
    ymgta=locsall.ynm_gt-y0gt;
    [thagta,rhoagta]=cart2pol(xmgta,ymgta);
    inrgt=rhoagta>R-dR&rhoagta<R+dR;
    coordinates_gta=struct('rho',rhoagta,'theta',thagta,'drho',rhoagta*0+5, 'dtheta',rhoagta*0+pi/64,'x',locsall.xnm_gt,'y',locsall.ynm_gt,'frame',locsall.frame);
    [numbercornerassigned_gta,~,direct_gta,timing_gt]=assigntocornersdirect(coordinates_gta,inrgt,8,[],timepoints,ax3b);
    [numgap_gta, numgapf_gta]=countspacing(coordinates_gta,inrgt,8,ax4);
%     [th_gt,rho_gt]=cart2pol(locsall.xnm_gt(:)-x0,locsall.ynm_gt(:)-y0);
%     [th_gtf,rho_gtf]=cart2pol(locs.xnm_gt(inr)-x0,locs.ynm_gt(inr)-y0);
%     tg=sort(th_gt);tg(end+1)=tg(1)+2*pi;
%     tgf=sort(th_gtf);tgf(end+1)=tgf(1)+2*pi;
%     numlocs=sum(diff(tg)>step*.7);
%     numlocsf=sum(diff(tgf)>step*.7);
%     out.numcornersfiltered_gt=mean([numbercornerassigned_gt,numgap_gt],'omitnan');
%     out.numcornersall_gt=mean([numbercornerassigned_gta,numgap_gta],'omitnan');

    out.numcornersfiltered_gt=direct_gt;
    out.numcornersall_gt=direct_gta;
    out.timing_gt=timing_gt;

    strt_gt=[', gt filtered (m/d): ' num2str(mean([numgap_gt, numbercornerassigned_gt],'omitnan')) ', ' num2str(direct_gt)];
    strt_gt=[strt_gt ', gt all (m/d): ' num2str(mean([numgap_gta, numbercornerassigned_gta],'omitnan')) ', ' num2str(direct_gta)];
    
else
    numlocsf=0;
    strt_gt='';
end
if obj.display
   ax= obj.setoutput('localizations');
%    bar(ax1,numberincorners);
%    xlabel(ax1,'number of locs per corner assigned');
%    tstr={['corners: ' int2str(numlocsf)  ', assigned: ' ,int2str(numbercornerassined)], ['from gaps: ' int2str(numfound) ', fractional: ' num2str(numfound2,3) ]};
%    title(ax1,tstr)
%    
%    ax3= obj.setoutput('gaps');
%    histogram(ax3,dth/step,0:.1:max(dth/step)+0.1);
%    xlabel(ax3,'length of gap (in units of 2pi/8)');
%     tstr={['corners: ' int2str(numlocsf)  ', from gaps: ' int2str(numfound) ', fractional: ' num2str(numfound2,3) ], ['assigned: ' ,int2str(numbercornerassined)]};
%  title(ax3,tstr)
 
%    ax2= obj.setoutput('localizations');
hold(ax,'off')
   plot(ax,locs.xnm(inr),locs.ynm(inr),'*')
      hold(ax,'on')
   if ~isempty(locs.xnm_gt)
        plot(ax,locsall.xnm_gt,locsall.ynm_gt,'ro',locs.xnm_gt(inr),locs.ynm_gt(inr),'r+')
   end
   circle(x0,y0,R,'Parent',ax);
   circle(x0,y0,R-dR,'Parent',ax,'LineStyle',':');     
   circle(x0,y0,R+dR,'Parent',ax,'LineStyle',':');
   
   circle(xfr,yfr,Rfr,'Parent',ax,'LineStyle','--');
   
   angles=mdt+[0:pi/4:pi];
   xt=[x0-R*1.5 x0+R*1.5]';
   xp=x0+(xt-x0)*cos(angles);yp=y0+(xt-x0).*sin(angles);
   plot(ax,xp,yp,'r')
   xp=x0+(xt-x0)*cos(angles+pi/8);yp=y0+(xt-x0).*sin(angles+pi/8);
   plot(ax,xp,yp,'k')
   str=['assigned: ' num2str(numbercornerassigned) ', direct: ' num2str(assigneddirect), ' , gap (i/f): ' num2str(numfound) ', ' num2str(numfound2,2) strt_gt];
   
   title(ax,str);
%    title(ax,['theta=' num2str(((mdt))/pi*180) ', pos=' num2str(x0-obj.site.pos(1)) ,',' num2str(y0-obj.site.pos(2))])
   axis(ax,'equal')
   xlim(ax,xt)
   ylim(ax,[y0-R*1.5 y0+R*1.5])
end

%    ax3= obj.setoutput('image');
% figure(188);
% subplot(3,1,1);bar(numberincorners);
% title(['theta=' num2str(((mdt))/pi*180) ', pos=' num2str(x0-obj.site.pos(1)) ,',' num2str(y0-obj.site.pos(2))])
% subplot(3,1,2);histogram(dthplot,-step/2:pi/64:step/2);
% 
% subplot(3,1,3)
% plot(locs.xnm(inr),locs.ynm(inr),'x')
% if ~isempty(locs.xnm_gt)
%     plot(locs.xnm(inr),locs.ynm(inr),'x',locsall.xnm_gt,locsall.ynm_gt,'ro',locs.xnm_gt(inr),locs.ynm_gt(inr),'r*')
% end
%    axis equal
%     title(['corners: ' int2str(numlocsf) ', found: ' int2str(numfound) ', fractional: ' num2str(numfound2,3), ', assigned: ' ,int2str(numbercornerassined)])
%     
% end
end

function [numbercornerassigned,mdt,direct]=assigntocorners(locs,inr,corners,ax,co)
if nargin <5 || isempty(co)
    co=0.5;
end
%assign to corneres
% find rotation
step=2*pi/corners;
% locptheta=locs.drho(inr)./locs.rho(inr);

mdt=cyclicaverage(locs.theta(inr),step,1./locs.dtheta(inr).^2);
if mdt>pi/16
    mdt=mdt-step;
end
dthplot=mod(locs.theta(inr)-mdt+step/2,step)-step/2;

%locptheta: increase for inefficient rotation etc:
spreadcorners=pi/32; %somehting like that
additionalerror=pi/64;
locpthetacorr=locs.dtheta(inr)+spreadcorners+additionalerror;
%somewhere: additional weighting with localiaztion precision: really bad
%ones: should not contribute
throt=locs.theta(inr)-mdt; %rotate with respect to template.
% cyclicaverage(throt,step,1./locptheta.^2)/pi*180 %funny: thats not zero.
% Evaluate to check goodness of rotation?
cornerpos=0:step:2*pi-step;
probc=zeros(8,length(throt));
probdirect=zeros(8,length(throt));
% cornerfills=zeros(8,1);
for k=1:length(throt)
    for c=1:8 %8 corners
        dc=mod(throt(k)-cornerpos(c)+pi,2*pi)-pi;
        probc(c,k)=exp(-dc^2/2/locpthetacorr(k)^2);
        if abs(dc)<step/2
            probdirect(k,c)=1;
        end
    end
%     probch=probc(:,k);
%     probch(probch<0.3)=0;
%     if sum(probch)<0.1
%         continue
%     end
%     probch=probch/sum(probch);
    
%     cornerfills=cornerfills+probch;
end
probc2=probc./sum(probc,1);
probc2=probc2/length(throt)*corners;
% probc2(probc2<co/2)=0;
numberincorners=sum(probc2,2);
numbercornerassigned=sum(numberincorners>co);
if numbercornerassigned==0
    numbercornerassigned=NaN;
end

direct=sum(any(probdirect,1));
if nargin>3 &&~isempty(ax)
bar(ax,numberincorners);
xlabel(ax,'number of locs per corner assigned');
tstr={['assigned: ' ,int2str(numbercornerassigned) ', direct: ' int2str(direct)]};
title(ax,tstr)
end
end


function [numbercornerassigned,mdt,direct,timing]=assigntocornersdirect(locs,inr,corners,ax,timepoints,axtp)
% if nargin <5 || isempty(co)
%     co=0.5;
% end
%assign to corneres
% find rotation
step=2*pi/corners;
% locptheta=locs.drho(inr)./locs.rho(inr);

mdt=cyclicaverage(locs.theta(inr),step,1./locs.dtheta(inr).^2);
if mdt>pi/16
    mdt=mdt-step;
end
frameh=locs.frame(inr);
throt=locs.theta(inr)-mdt; %rotate with respect to template.
throt=mod(throt-step/2,2*pi);
cornerposd=0:step:2*pi;
cornerposdf=0:step/8:2*pi;
h=histcounts(throt,cornerposd);
hf=histcounts(throt,cornerposdf);
direct=sum(h>=1);


numbercornerassigned=direct;
if nargin>3 &&~isempty(ax)
     hold(ax,'off')
    bar(ax,(cornerposd(1:end-1)+step/2)/step,h)
     hold(ax,'on')
     bar(ax,(cornerposdf(1:end-1)+step/2/8)/step,hf)
%      histogram(ax,throt/step,64)
% bar(ax,numberincorners);
xlabel(ax,'number of locs per corner assigned');
tstr={['assigned direct: ' int2str(direct)]};
title(ax,tstr)
end


%temporal dependence
if nargin>4 && ~isempty(timepoints)
    nstart=zeros(length(timepoints),1);
    nend=zeros(length(timepoints),1);
    for k=1:length(timepoints)
        indin=frameh<timepoints(k);
        nstart(k)=sum(histcounts(throt(indin),cornerposd)>=1);
        nend(k)=sum(histcounts(throt(~indin),cornerposd)>=1);
    end
    if nargin>5 &&~isempty(axtp)
        plot(axtp,timepoints,nstart,'-*',timepoints,nend,'-o')
    end
    timing.timepoints=timepoints;
    timing.nstart=nstart;
    timing.nend=nend;
else
    timing=[];
end

end

function [numfound, numfound2]=countspacing(locs,inr,corners,ax)
step=2*pi/corners;
th=locs.theta(inr);
locptheta=locs.dtheta(inr);
%number of locs from spacing
if isempty(locptheta)
    numfound=0;numfound2=0;
    return
end

[ths,inds]=sort(th);
ths(end+1)=ths(1)+2*pi;
dth=diff(ths);
lp=locptheta(inds);
lp(end+1)=lp(1);
dstep=sqrt(lp(1:end-1).^2+lp(2:end).^2); %uncertaintty
mindistance=1;%-0.1*length(th)/8;

fullspace=ceil(dth-mindistance); fullspace( fullspace<0)=0;
numfound=8-sum(fullspace);

% try again with non-integer values
fs2=dth; fs2=fs2(fs2>mindistance);
fs2=fs2-.25;
numfound2=8-sum(fs2);

if nargin>3 && ~isempty(ax)
 histogram(ax,dth/step,0:.1:max(dth/step)+0.1);
   xlabel(ax,'length of gap (in units of 2pi/8)');
    tstr={[ 'from gaps: ' int2str(numfound) ', fractional: ' num2str(numfound2,3) ]};
 title(ax,tstr)
%  hold(ax,'on')
 
end
end