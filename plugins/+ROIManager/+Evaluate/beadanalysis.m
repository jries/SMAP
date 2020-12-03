classdef beadanalysis<interfaces.SEEvaluationProcessor
%     Calculates the line-profile along user-defined direction and fits it
%     with a Gaussian or double-Gaussian model
    properties
%         boundary
    end
    methods
        function obj=beadanalysis(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj, p)
        
            locs=obj.getLocs({'xnm','ynm','frame','phot'},'layer',1,'size',p.se_siteroi/2);  
            ax=obj.setoutput('intensity');
            plot(ax,locs.frame,locs.phot)
            
   
            out.frame=locs.frame;
            out.phot=locs.phot;
            out.noise=std(diff(diff(locs.phot)))/2;
        end
     
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end

end



function pard=guidef(obj)
% pard.fitmodel.object=struct('Style','popupmenu','String',{{'Gaussian','Step','Disk','Ring','Double Gaussian'}});
% pard.fitmodel.position=[1,1];
% pard.fitmodel.Width=2;
% 
% pard.restrictsigma.object=struct('String','sigma=<locp>','Style','checkbox');
% pard.restrictsigma.position=[3,1];
% pard.restrictsigma.Width=4;
%  p(1).value=0; p(1).on={}; p(1).off={'binwidth'};
% p(2).value=1; p(2).on={'binwidth'}; p(2).off={};
% pard.setbinwidth.object=struct('String','set binwidth (nm):','Style','checkbox','Value',1,'Callback',{{@obj.switchvisible,p}});
% pard.setbinwidth.position=[2,1];
% pard.setbinwidth.Width=3;
% pard.setbinwidth.Tooltip='Set bin width for profile. If not set, use pixel size from ROI manager';
% pard.binwidth.object=struct('String','2','Style','edit');
% pard.binwidth.position=[2,3.5];
% pard.binwidth.Width=0.5;
% pard.binwidth.TooltipString='Binwidth for profiles. If not checked, use pixel size of reconstruction';
% pard.setbinwidth.TooltipString=pard.binwidth.TooltipString;
% 
% % pard.dxt.Width=3;
pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi','se_sitepixelsize'};
pard.plugininfo.type='ROI_Evaluate';
pard.plugininfo.description='beads';
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