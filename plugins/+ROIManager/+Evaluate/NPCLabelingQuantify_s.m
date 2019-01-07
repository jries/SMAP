classdef NPCLabelingQuantify_s<interfaces.SEEvaluationProcessor
    properties
        savedevals
    end
    methods
        function obj=NPCLabelingQuantify_s(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function makeGui(obj)
            makeGui@interfaces.SEEvaluationProcessor(obj);
        end
        function out=run(obj,p)
            out=runintern(obj,p);         
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

pard.minlocst.object=struct('Style','text','String','min locs/corner');
pard.minlocst.position=[2,1];
pard.minlocst.Width=3;

pard.minlocs.object=struct('Style','edit','String','1');
pard.minlocs.position=[2,4];
pard.minlocs.Width=1;

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
out.numlocsf=sum(inr);
out.numlocsR=sum(rhoa<R+dR);
out.coordinates=struct('rho',rhoa,'theta',tha,'drho',locs.locprecnm, 'dtheta',locs.locprecnm./rhoa,'x',locs.xnm,'y',locs.ynm,'frame',locs.frame);


if obj.display
   axdat=obj.setoutput('Data');
   tabdat=uitabgroup(axdat.Parent);
   ax1=axes(uitab(tabdat,'Title','assigned'));
   ax1b=axes(uitab(tabdat,'Title','a_t'));
   if ~isempty(locs.xnm_gt)
        axgt=obj.setoutput('Ground_truth');
       tabdat=uitabgroup(axgt.Parent);
       ax3=axes(uitab(tabdat,'Title','assigned'));
       ax3b=axes(uitab(tabdat,'Title','a_t'));
   end
else
    ax1=[]; ax3=[]; ax3b=[]; ax1b=[];
end
  
locsall=obj.locData.getloc({'frame','filenumber'},'layer',1,'removefilter','filenumber','position','all');

frames=locsall.frame(locsall.filenumber==obj.site.info.filenumber);
frames0=frames(frames>0);
numpoints=50;
qq=linspace(0,1,numpoints+1);
timepoints=myquantile(frames0,qq);  % now edges are defined: think about time windows

[numbercornerassigned,mdt,assigneddirect,timing]=assigntocornersdirect(out.coordinates,inr,8,ax1,timepoints,ax1b,p.minlocs);
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
    [numbercornerassigned_gt,~,direct_gt,timing_gt]=assigntocornersdirect(coordinates_gt,inr,8,ax3,timepoints);

    locsall=obj.getLocs({'xnm','ynm','xnm_gt','ynm_gt','locprecnm','frame'},'size',p.se_siteroi);
    xmgta=locsall.xnm_gt-x0gt;
    ymgta=locsall.ynm_gt-y0gt;
    [thagta,rhoagta]=cart2pol(xmgta,ymgta);
    inrgt=rhoagta>R-dR&rhoagta<R+dR;
    coordinates_gta=struct('rho',rhoagta,'theta',thagta,'drho',rhoagta*0+5, 'dtheta',rhoagta*0+pi/64,'x',locsall.xnm_gt,'y',locsall.ynm_gt,'frame',locsall.frame);
    [numbercornerassigned_gta,~,direct_gta]=assigntocornersdirect(coordinates_gta,inrgt,8,[],timepoints,ax3b);

    out.numcornersfiltered_gt=direct_gt;
    out.numcornersall_gt=direct_gta;
    out.timing_gt=timing_gt;

    strt_gt=[', gt filtered (m/d): ' num2str(mean([numbercornerassigned_gt],'omitnan')) ', ' num2str(direct_gt)];
    strt_gt=[strt_gt ', gt all (m/d): ' num2str(mean([ numbercornerassigned_gta],'omitnan')) ', ' num2str(direct_gta)];
    
else
    strt_gt='';
end
if obj.display
   ax= obj.setoutput('localizations');
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
   str=['assigned: ' num2str(numbercornerassigned) ', direct: ' num2str(assigneddirect)];
   
   title(ax,str);
   axis(ax,'equal')
   xlim(ax,xt)
   ylim(ax,[y0-R*1.5 y0+R*1.5])
end
end


function [numbercornerassigned,mdt,direct,timing]=assigntocornersdirect(locs,inr,corners,ax,timepoints,axtp,minlocs)
if nargin<7
    minlocs=1;
end

%assign to corneres
% find rotation
step=2*pi/corners;
mdtc=cyclicaverage(locs.theta(inr),step,1./locs.dtheta(inr).^2);
if mdtc>pi/16
    mdtc=mdtc-step;
end

% with fitting
thh=double(locs.theta(inr));

% ft=fittype('mod(x-a,pi/4)');
% fitp=fit(thh,0*thh+pi/8,ft,'Weights',double(1./locs.dtheta(inr)),'StartPoint',mdtc+pi/8);

%weighted fit
fitfun=@(x,xdata)mod(x-xdata,pi/4);
y=pi/8;
weights=double(1./locs.dtheta(inr));
costfun=@(A) weights.*(fitfun(A,thh)-y);
options=optimoptions('lsqnonlin','Display','off');
fitpl=lsqnonlin(costfun,double(mdtc)+pi/8,[],[],options);
mdt=fitpl-pi/8;
% mdt=fitp.a-pi/8;


frameh=locs.frame(inr);
throt=locs.theta(inr)-mdt; %rotate with respect to template.
throt=mod(throt-step/2,2*pi);
cornerposd=0:step:2*pi;
cornerposdf=0:step/8:2*pi;
h=histcounts(throt,cornerposd);
hf=histcounts(throt,cornerposdf);
direct=sum(h>=minlocs);

numbercornerassigned=direct;
if nargin>3 &&~isempty(ax)
     hold(ax,'off')
    bar(ax,(cornerposd(1:end-1)+step/2)/step,h)
     hold(ax,'on')
     bar(ax,(cornerposdf(1:end-1)+step/2/8)/step,hf)
    xlabel(ax,'number of locs per corner assigned');
    tstr={['assigned direct: ' int2str(direct)]};
    title(ax,tstr)
end


%temporal dependence
if nargin>4 && ~isempty(timepoints)
    nstart=zeros(length(timepoints),1);
    nend=zeros(length(timepoints),1);
    for k=1:length(timepoints)
        indin=frameh<=timepoints(k);
        nstart(k)=sum(histcounts(throt(indin),cornerposd)>=1);
        nend(k)=sum(histcounts(throt(~indin),cornerposd)>=1);
    end
    if nargin>5 &&~isempty(axtp)
        plot(axtp,timepoints,nstart,'-*',timepoints,nend,'-o')
    end
    timing.timepoints=timepoints;
    timing.nstart=nstart;
    timing.nend=nend;
    
    % new analysis: divide into chunks
    ltp=(length(timepoints))-1;
    nchunks=zeros(ltp,ltp)+NaN;nchunksn=nchunks;
    dind=nchunks;
    for d=1:ltp
        for k=1:ltp
            k2=mod(k+d,ltp+1); 
            if k2==0 
                k2=ltp+1; 
            end
            if k2>k
                indin=frameh>=timepoints(k) & frameh<timepoints(k2);
            else
                indin=frameh>=timepoints(k) | frameh<timepoints(k2+1);
            end
            
            if sum(indin)>0
                nchunks(d,k)=sum(histcounts(throt(indin),cornerposd)>=1);
%                 nchunksn(d,k)=sum(histcounts(throt(~indin),cornerposd)>=1);
            else
                 nchunks(d,k)=0;
            end
             dind(d,k)=d;
        end    
    end
    timing.nchunks=nchunks;
    timing.nchunksn=nchunksn;
    timing.dind=dind;
else
    timing=[];
end

end

