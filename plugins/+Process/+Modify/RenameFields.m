classdef RenameFields<interfaces.DialogProcessor
    methods
        function obj=RenameFields(varargin)     
            obj@interfaces.DialogProcessor(varargin{:});
        end
        
        function out=run(obj,p)
            obj.setPar('undoModule','RemoveFields');
            notify(obj.P,'backup4undo');
            oldfield=p.fieldselect.selection;
            if isempty(p.newfield)
                obj.locData.loc=myrmfield(obj.locData.loc,oldfield);
            else
            
            %if new exists: rename it by adding "_old"
            
            if isfield(obj.locData.loc,p.newfield)
                obj.locData.loc.([p.newfield '_old'])=obj.locData.loc.(p.newfield);
            end
            
            %clear filters of old
            for k=1:obj.getPar('numberOfLayers')
                obj.locData.layer(k).filter=myrmfield(obj.locData.layer(k).filter,oldfield);
            end
            
           %rename old to new 
            obj.locData.loc.(p.newfield)=obj.locData.loc.(oldfield);
            if ~p.keepold
                obj.locData.loc=myrmfield(obj.locData.loc,oldfield);
            end
            end
            obj.locData.regroup;
            obj.setPar('locFields',fieldnames(obj.locData.loc))  
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

pard.textb.object=struct('String','old field','Style','text');
pard.textb.position=[2,1];
pard.textb.object.TooltipString='select field which you want to rename';


pard.fieldselect.object=struct('Style','popupmenu','String','field');
pard.fieldselect.position=[2,2];
pard.fieldselect.object.TooltipString='select field which you want to rename';
pard.textb.object.TooltipString=pard.fieldselect.object.TooltipString;

pard.textc.object=struct('String','new field name','Style','text');
pard.textc.position=[3,1];

pard.newfield.object=struct('Style','edit','String','field');
pard.newfield.position=[3,2];
pard.newfield.object.TooltipString='new field name. Leave empty to delete field';
pard.textc.object.TooltipString=pard.newfield.object.TooltipString;


pard.keepold.object=struct('Style','checkbox','String','keep old field','Value',1);
pard.keepold.position=[4,1];
pard.keepold.object.TooltipString='If checked: keep old field (duplicate), if unchecked: delete old field';

pard.syncParameters={{'locFields','fieldselect',{'String'}}};
pard.plugininfo.type='ProcessorPlugin';
end