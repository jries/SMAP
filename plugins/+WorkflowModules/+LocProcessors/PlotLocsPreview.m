classdef PlotLocsPreview<interfaces.WorkflowModule;
    properties
%          pixelsize      
    end
    methods
       function obj=PlotLocsPreview(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=1; 
%              obj.setInputChannels(1,'frame');
       end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
        end
        function prerun(obj,p)
            
        end
        function pard=guidef(obj)
            pard.plugininfo.type='WorkflowModule'; 
            pard.plugininfo.description='Plots localizations during preview.';
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
                    for k=1:length(locs.xpix)
                        pos=[locs.xpix(k)-dn locs.ypix(k)-dn locs.xpix(k)+dn locs.ypix(k)+dn ];
                        col=[1 0. 1];
                        plotrect(ax,pos,col);
%                         plotrect(ax,pos,{'Color',col,'LineWidth',3});
                    end
      
                end
            end 


        end

    end
end
