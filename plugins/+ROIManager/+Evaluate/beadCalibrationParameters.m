classdef beadCalibrationParameters<interfaces.SEEvaluationProcessor
    properties
        
    end
    methods
        function obj=beadCalibrationParameters(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
                obj.inputParameters={'se_cellfov','se_sitefov','se_siteroi'};
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
function out=runintern(obj,p)

ax=obj.setoutput('sx');
hold(ax, 'off')
ax2=obj.setoutput('sy');
hold(ax2, 'off')
axp=obj.setoutput('phot');
hold(axp, 'off')
sites=obj.locData.SE.sites;
for s=length(sites):-1:1
    if sites(s).annotation.use
        posh=horzcat(sites(s).pos(1:2),p.se_siteroi);
        locs=obj.locData.getloc({'xnm','ynm','PSFxnm','PSFynm','frame','phot'},'Position',posh,'layer',1,'grouping','ungrouped');
        
        plot(ax, locs.frame,locs.PSFxnm,'-')
        hold(ax, 'on')
        plot(ax2, locs.frame,locs.PSFynm,'-')
        hold(ax2, 'on')      
        plot(axp, locs.frame,locs.phot,'-')
        hold(axp, 'on')
       
    end
    
end

        posh=horzcat(obj.site.pos(1:2),p.se_siteroi);
        locs=obj.locData.getloc({'xnm','ynm','PSFxnm','PSFynm','frame','phot'},'Position',posh,'layer',1,'grouping','ungrouped');
        
        plot(ax, locs.frame,locs.PSFxnm,'k-')
        hold on
        plot(ax2, locs.frame,locs.PSFynm,'k-')
        hold on        
        plot(axp, locs.frame,locs.phot,'k-')
        hold on
        out=[];
end

function pard=guidef
% pard.centermode.object=struct('Style','popupmenu','String',{{'median','mean','mask','fitring'}});
% pard.centermode.position=[1,1];
% pard.centermode.Width=2;
% 
% pard.iterationst.object=struct('Style','text','String','iterations: ');
% pard.iterationst.position=[2,1];
% pard.iterationst.Width=1;
% 
% pard.iterations.object=struct('Style','edit','String',3);
% pard.iterations.position=[2,2];
% pard.iterations.Width=1;
% 
% pard.cxy.object=struct('Style','checkbox','String','correct x,y');
% pard.cxy.position=[1,3];
% pard.cxy.Width=2;
% 
% pard.cz.object=struct('Style','checkbox','String','correct z');
% pard.cz.position=[2,3];
% pard.cz.Width=2;
% 
 pard.plugininfo.type='ROI_Evaluate';
% pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi','layer1_','layer2_','se_sitepixelsize'};
end

