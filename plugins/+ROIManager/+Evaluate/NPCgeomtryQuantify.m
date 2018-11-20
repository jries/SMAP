classdef NPCgeomtryQuantify<interfaces.SEEvaluationProcessor
    properties
        savedevals
    end
    methods
        function obj=NPCgeomtryQuantify(varargin)        
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

pard.R.object=struct('Style','edit','String','50');
pard.R.position=[1,2];
pard.R.Width=1;

pard.dRt.object=struct('Style','text','String','dR:');
pard.dRt.position=[1,3];
pard.dRt.Width=1;

pard.dR.object=struct('Style','edit','String','20');
pard.dR.position=[1,4];
pard.dR.Width=1;

pard.plugininfo.type='ROI_Evaluate';
pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi','layer1_','layer2_','se_sitepixelsize'};
end


function out=runintern(obj,p)
R=p.R;
dR=p.dR;

locs=obj.getLocs({'xnm','ynm','znm','locprecnm','locprecznm','frame'},'layer',1,'size',p.se_siteroi(1)/2);
if isempty(locs.xnm)
    out=[];
    return
end


[x0,y0]=fitposring(locs.xnm,locs.ynm,R);
xm=locs.xnm-x0;
ym=locs.ynm-y0;

if ~isempty(locs.znm) % evaluate z distance
    dz=8;
    z0=median(locs.znm);
    z=-100:dz:100+z0;
    hz=hist(locs.znm-obj.site.pos(3),z);
    % hz=hz-mean(hz);
    ac=myxcorr(hz,hz);
    % if obj.display
    if obj.display
    ax1=obj.setoutput('profile');
       fitresult=createFit(z, hz,ax1);
    title(ax1,fitresult.d)
     ax2=obj.setoutput('correlation');
    plot(ax2,z,ac)
    else
    fitresult=createFit(z, hz,[]);
    end
    
    % plot(ax1,z,hz)
  


    % end
    out.ac=ac;
    out.dz=dz;
    out.hist=hz;
    out.z=z;
    out.Gaussfit=copyfields([],fitresult,{'a1','a2','b','c','d'});
end

%radial analysis
thetan=-pi:pi/128:pi;
[theta,rhoa]=cart2pol(xm,ym);
histtheta=hist(theta,thetan);
% histtheta12=hist(theta+pi,thetan+pi);
% tac=myxcorr(histtheta11,histtheta11);
% tac1=tac+myxcorr(histtheta12,histtheta12);
tac1=xcorrangle(histtheta-mean(histtheta));
if obj.display
   ax1=obj.setoutput('corrtheta');
   plot(ax1,thetan-thetan(1),tac1);
end

out.actheta=tac1;
out.thetan=thetan-thetan(1);
end


function [fitresult, gof] = createFit(z, hz,ax)

[xData, yData] = prepareCurveData( z, hz );

% Set up fittype and options.
ft = fittype( 'a1*exp(-((x-b)/c)^2/2) + a2*exp(-((x-b+d)/c)^2/2)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 -Inf 0 -Inf];
mv= sum(hz.*z)/sum(hz);
dh=30;
sp=[max(hz) max(hz) mv+dh 8 dh*2 ];
opts.StartPoint = sp;
fstart=ft(sp(1),sp(2),sp(3),sp(4),sp(5),xData);
% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
if ~isempty(ax)
axes(ax)
h = plot(xData,fitresult(xData),'r', xData, yData,'b-',xData,fstart,'g');


% legend( h, 'hz vs. z', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
xlabel z
ylabel hz
grid on
end

end
