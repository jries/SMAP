classdef SortCells<interfaces.DialogProcessor&interfaces.SEProcessor
    methods
        function obj=SortCells(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'se_viewer'};
        end
        
        function out=run(obj,p)  
            out=[];
%             sites=obj.locData.SE.sites;
            cells=obj.locData.SE.cells;
            nl=zeros(length(cells),1);
            fn=zeros(length(cells),1);
            cn=zeros(length(cells),1);
            for k=1:length(cells)
                nl(k)= cells(k).image.layers(1).images.rawimage.numberOfLocs;
                fn(k)=cells(k).info.filenumber;
                cn(k)=cells(k).ID;
            end
%             nosites=false(length(obj.SE.cells),1);
%             info=[sites(:).info];
%             sitecells=[info.cell];
%             for k=1:length(cells)
%                 if ~any(cells(k).ID==sitecells)
%                     obj.locData.SE.removeCell(cells(k).ID);
%                 end
%             end
            switch p.sortselection.Value
                case 1 %nl
                    [~,indsort]=sort(nl);
                case 2
                    [~,indsort]=sort(cn);
                case 3
                    [~,indsort]=sort(fn);
            end
            
            
            cells=cells(indsort);
            obj.locData.SE.cells=cells;
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

pard.sortselection.object=struct('String',{{'number of localizations','Cell number' 'File number'}},'Style','popupmenu');
pard.sortselection.position=[1,1];
pard.sortselection.Width=2;

pard.plugininfo.type='ROI_Analyze';


end