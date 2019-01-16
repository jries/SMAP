classdef lineprofile<interfaces.SEEvaluationProcessor
    properties
        boundary
    end
    methods
        function obj=lineprofile(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj, inp)
            % Import the kimograph
%             site = obj(1).locData
%             runNPC3Dfitting(obj,inp)
            layers=find(inp.sr_layerson);
            ax=obj.setoutput('profile');
            for k=1:length(layers)
                locs=obj.getLocs({'xnmrot','ynmrot','znm'},'layer',layers(k),'size',inp.se_siteroi/2);  
                nbins=-200:5:200;
                hc=histcounts(locs.ynmrot,nbins);
                mp=median(locs.ynmrot);
                d=std(locs.ynmrot);
%                 d=30;
                amp=max(hc);
                ft=fittype( @(a1,b1,s,a2,b2,x) a1*exp(-(x-b1).^2/2/s^2)+a2*exp(-(x-b2).^2/2/s^2));
%                 ft=fittype( @(a1,b1,s,a2,b2,c,x) a1*normcdf((x-b1)/s)+a2*normcdf((x-b2)/s)+c);
%                 hc=cumsum(hc);
                fitp=fit(nbins(1:end-1)',hc',ft,'StartPoint',[amp mp-d d amp mp+d]);
                plot(ax,nbins(1:end-1),hc,nbins,fitp(nbins));
                dist=abs(fitp.b2-fitp.b1);
                title(ax,dist)

            end
            out.distance=dist;
            out.std=fitp.s;
        end
     
        function pard=guidef(obj)
            pard=guidef;
        end
    end

end



function pard=guidef
pard.fitmodel.object=struct('Style','popupmenu','String',{{'Gaussian','Double Gaussian'}});
pard.fitmodel.position=[1,1];
pard.fitmodel.Width=2;


% pard.dxt.Width=3;
pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi'};
pard.plugininfo.type='ROI_Evaluate';
end

function M = calMeasurement(x,y,qx,qy,Size)
    ref = interp1(x,y,qx);
    rightIdx = qy > ref;
    leftIdx = ~rightIdx;
    leftA = sum(y);
    rightA = Size(1)*Size(2)-sum(y);
    leftD = sum(leftIdx)/leftA;
    rightD = sum(rightIdx)/rightA;
    M = 0-((leftD-1)^2 + (rightD-0)^2)^(1/2);
end