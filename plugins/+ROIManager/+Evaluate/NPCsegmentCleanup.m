classdef NPCsegmentCleanup<interfaces.SEEvaluationProcessor
    properties
        savedevals
    end
    methods
        function obj=NPCsegmentCleanup(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function makeGui(obj)
            makeGui@interfaces.SEEvaluationProcessor(obj);
%             obj.guihandles.saveimagesb.Callback={@saveimagesb_callback,obj};
        end
        function out=run(obj,p)
            try
            out=runintern(obj,p);
            catch err
                err
                out=[];
            end
         
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

pard.R.object=struct('Style','edit','String','55');
pard.R.position=[1,2];
pard.R.Width=1;

pard.dRt.object=struct('Style','text','String','dR:');
pard.dRt.position=[1,3];
pard.dRt.Width=1;

pard.dR.object=struct('Style','edit','String','15');
pard.dR.position=[1,4];
pard.dR.Width=1;

% pard.remove.object=struct('Style','checkbox','String','remove bad');
% pard.remove.position=[6,1];
% pard.remove.Width=2;

pard.center.object=struct('Style','checkbox','String','center','Value',1);
pard.center.position=[2,1];
pard.center.Width=1;

pard.insidet.object=struct('Style','text','String','in/ring <');
pard.insidet.position=[3,1];
pard.insidet.Width=1.5;
pard.inside.object=struct('Style','edit','String','0.3');
pard.inside.position=[3,2];
pard.inside.Width=.5;

pard.outsidet.object=struct('Style','text','String','out/ring <');
pard.outsidet.position=[3,3];
pard.outsidet.Width=1.5;
pard.outside.object=struct('Style','edit','String','.75');
pard.outside.position=[3,4];
pard.outside.Width=.5;

pard.radiusranget.object=struct('Style','text','String','Radius min max');
pard.radiusranget.position=[4,1];
pard.radiusranget.Width=2;
pard.radiusrange.object=struct('Style','edit','String','45 65');
pard.radiusrange.position=[4,3];
pard.radiusrange.Width=1;

pard.minsizet.object=struct('Style','text','String','min size');
pard.minsizet.position=[5,1];
pard.minsizet.Width=1.5;
pard.minsize.object=struct('Style','edit','String','25');
pard.minsize.position=[5,2.5];
pard.minsize.Width=0.5;

pard.minlocst.object=struct('Style','text','String','min localizations');
pard.minlocst.position=[5,3];
pard.minlocst.Width=1.5;
pard.minlocs.object=struct('Style','edit','String','10');
pard.minlocs.position=[5,4.5];
pard.minlocs.Width=.5;

pard.maxPSFt.object=struct('Style','text','String','max average PSF');
pard.maxPSFt.position=[6,1];
pard.maxPSFt.Width=1.5;
pard.maxPSF.object=struct('Style','edit','String','150');
pard.maxPSF.position=[6,2.5];
pard.maxPSF.Width=.5;


pard.plugininfo.type='ROI_Evaluate';
pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi','layer1_','layer2_','se_sitepixelsize'};
end


function out=runintern(obj,p)
R=p.R;
dR=p.dR;
    out=[];



%2D fit


if p.center %directly center to do furhter analysis on centered pore
    extraspace=1.2;
    locs=obj.getLocs({'xnm','ynm'},'layer',1,'size',p.se_siteroi(1)/2*extraspace);
    [x0,y0,R0,resnorm]=fitposring(locs.xnm,locs.ynm,R);
    obj.site.pos(1:2)=[x0, y0];

%     obj.redraw;
end

locs=obj.getLocs({'xnm','ynm','PSFxnm'},'layer',1,'size',p.se_siteroi(1)/2);
if isempty(locs.xnm)
    return
end

    [x0,y0,R0,resnorm]=fitposring(locs.xnm,locs.ynm,R);

xm=locs.xnm-x0;
ym=locs.ynm-y0;
[x0r,y0r,R0r,resnormR]=fitposring(xm,ym);



if obj.display
    ax1=obj.setoutput('ringfit');

    plot(ax1,xm-x0r,ym-y0r,'.')
    hold(ax1,'on')
    circle(0,0,R0r,'Parent',ax1);
    circle(-x0r,-y0r,R0,'Parent',ax1,'LineStyle','--');
    hold(ax1,'off')
    title(ax1,R0r)
end

% fitted radius
goodradius=R0r>p.radiusrange(1) & R0r<p.radiusrange(2);
% in ring vs background
insidei=(xm.^2+ym.^2<(R-dR)^2);
inringi=(xm.^2+ym.^2>(R-dR)^2 & xm.^2+ym.^2<(R+dR)^2);
outsidei=(xm.^2+ym.^2>(R+dR)^2);
inside=sum(insidei);
inring=sum(inringi);
outside=sum(outsidei);
goodin=inside/inring<p.inside;
goodout=outside/inring<p.outside;

% minimum size 
sm=sqrt(eig(cov(xm,ym)));
smm=min(sm);
smav=sqrt(prod(sm));

goodsize=smav>p.minsize;
goodlocs=length(xm)>p.minlocs;
mpsf=mean(locs.PSFxnm);

goodpsf=mpsf<=p.maxPSF;
usethis=goodradius&goodout&goodin&goodsize&goodlocs & goodpsf;
    savefield='list4';
obj.site.annotation.(savefield).value=usethis+1;
obj.site.annotation.use=usethis;

out.R0=R0r;out.inside=inside;out.outside=outside;out.inring=inring;
out.sizeav=smav;
out.locs=length(xm);

if obj.display
    textstyle={'\rm','\bf'};
    textstart={'\rm bad: ', '\bf good: '};

    ax1=obj.setoutput('ringfit');

    plot(ax1,xm(inringi)-x0r,ym(inringi)-y0r,'.')
    hold(ax1,'on')
    plot(ax1,xm(~inringi)-x0r,ym(~inringi)-y0r,'r.')
    circle(0,0,R0r,'Parent',ax1,'EdgeColor','r');
    circle(-x0r,-y0r,R0,'Parent',ax1);
    circle(-x0r,-y0r,R0+dR,'Parent',ax1,'LineStyle','--');
    circle(-x0r,-y0r,R0-dR,'Parent',ax1,'LineStyle','--');
    hold(ax1,'off')

    ff='%2.1f';
    ff2='%1.2f';
    titletxt=[textstart{usethis+1} '\bf R: ' textstyle{goodradius+1} num2str(R0r,ff) ...
        '\bf, in: ' textstyle{goodin+1} num2str(inside/inring,ff2)...
        '\bf, out: ' textstyle{goodout+1} num2str(outside/inring,ff2)...
        '\bf, size: ' textstyle{goodsize+1} num2str(smav,ff) ...
        '\bf, locs: ' textstyle{goodlocs+1} num2str(length(xm),'%2i') ...
        '\bf, psf: ' textstyle{goodpsf+1} num2str(mpsf,'%3.0f')];
    title(ax1,titletxt)
end



% minimum localizations
end

