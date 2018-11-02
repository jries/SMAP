classdef RemoveEmptyCells<interfaces.DialogProcessor&interfaces.SEProcessor
    methods
        function obj=RemoveEmptyCells(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'se_viewer'};
        end
        
        function out=run(obj,p)  
            out=[];
            sites=obj.locData.SE.sites;
            cells=obj.locData.SE.cells;
%             nosites=false(length(obj.SE.cells),1);
            info=[sites(:).info];
            sitecells=[info.cell];
            for k=1:length(cells)
                if ~any(cells(k).ID==sitecells)
                    obj.locData.SE.removeCell(cells(k).ID);
                end
            end
            viewer=p.se_viewer;
            viewer.updateCelllist;
            if isempty(obj.locData.SE.indexFromID(cells,obj.locData.SE.currentcell.ID))

            obj.locData.SE.currentcell=obj.locData.SE.cells(1);
            end
                
           
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end




function pard=guidef
pard.plugininfo.type='ROI_Analyze';


end