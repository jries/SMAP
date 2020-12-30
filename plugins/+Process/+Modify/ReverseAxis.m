classdef ReverseAxis<interfaces.DialogProcessor&interfaces.SEProcessor
    % Reverse the specified axis across the entire file
    methods
        function obj=ReverseAxis(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
%             obj.showresults=false;
        end
        
        function out=run(obj,p)
            out=[];
            obj.setPar('undoModule','ReverseAxis');
            notify(obj.P,'backup4undo');
            
            %% Get info from locs
            locs = obj.locData.loc;
            selectedAxis = p.axis.selection;
            minOneAxis = min(locs.(selectedAxis));
            locs.(selectedAxis) = -locs.(selectedAxis);
            newMinOneAxis = min(locs.(selectedAxis));
            dMin = minOneAxis-newMinOneAxis;
            locs.(selectedAxis) = locs.(selectedAxis)+dMin;
            
            %% Apply changes to sites
            if obj.SE.numberOfSites>0
                sites = obj.SE.sites;
                switch (selectedAxis)
                    case 'xnm'
                        onePos = 1;
                    case 'ynm'
                        onePos = 2;
                end
                for k = 1:obj.SE.numberOfSites
                    sites(k).pos(onePos) = -sites(k).pos(onePos)+dMin;
                end
            end
            obj.locData.loc.(selectedAxis) = locs.(selectedAxis);           
            obj.locData.regroup;
            obj.locData.filter;
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end


function pard=guidef

pard.t.object=struct('String','Axis:','Style','text');
pard.t.position=[2 1];
pard.t.Width=1;

pard.axis.object=struct('Style','popupmenu','String',{{'xnm' 'ynm'}}, 'Value',1);
pard.axis.position=[2 2];
pard.axis.Width=1;

pard.plugininfo.type='ProcessorPlugin';

end