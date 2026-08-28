classdef FilterFields2<interfaces.DialogProcessor
%     Renames localization fields
    methods
        function obj=FilterFields2(varargin)     
            obj@interfaces.DialogProcessor(varargin{:});
        end
        
        function out=run(obj,p)
            out=[];
            % p.layers
            for k=1:length(p.layers)
                sf={p.fieldselect,p.filterminmax(1),p.filterminmax(2),1,1};
                obj.setPar('selectedField',sf,'layer',p.layers(k))
            end

            obj.locData.regroup;  
        end
        function pard=guidef(obj)
            pard=guidef;
        end

        function updateGui(obj,event,data)
            if ~isempty(obj.locData.loc)
            fn=fieldnames(obj.locData.loc);
            obj.guihandles.fieldselect.String=fn;
            end
        end       
    end
end




function pard=guidef

pard.textb.object=struct('String','Field name:','Style','text');
pard.textb.position=[1,1];

pard.fieldselect.object=struct('Style','edit','String','');
pard.fieldselect.position=[1,2];


pard.textc.object=struct('String','min,max','Style','text');
pard.textc.position=[1,3];

pard.filterminmax.object=struct('Style','edit','String','0,1');
pard.filterminmax.position=[1,4];


pard.textd.object=struct('String','Layers:','Style','text');
pard.textd.position=[2,1];

pard.layers.object=struct('Style','edit','String','1');
pard.layers.position=[2,2];


pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description='Apply localization filter';
end