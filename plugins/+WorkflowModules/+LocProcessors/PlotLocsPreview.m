classdef PlotLocsPreview<interfaces.WorkflowModule
%    Plots localizations of a single frame for optimization of fitting
%    parameters
    properties     
    end
    methods
       function obj=PlotLocsPreview(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=1; 
       end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
        end
        function prerun(obj,p)
            
        end
        function pard=guidef(obj)
            pard.plugininfo.type='WorkflowModule'; 
            pard.plugininfo.description='Plots localizations of a single frame for optimization of fitting parameters';
        end
        function output=run(obj,data,p)
            output=[];
            if obj.getPar('loc_preview') 
                locs=data.data;
                if ~isempty(locs)&&~isempty(locs.xpix)
                    ax=findobj(obj.getPar('loc_outputfig').Children,'Type','Axes');
                    ax.NextPlot='add';
                    plot(locs.xpix,locs.ypix,'k.','Parent',ax,'MarkerSize',4)
                    dn=ceil((obj.getPar('loc_ROIsize')-1)/2);
                    if isempty(dn)
                        dn=3;
                    end
                    for k=1:length(locs.xpix)
                        pos=[locs.xpix(k)-dn locs.ypix(k)-dn locs.xpix(k)+dn locs.ypix(k)+dn ];
                        col=[1 0. 1];
                        plotrect(ax,pos,col);
                    end
      
                end
            end 


        end

    end
end
