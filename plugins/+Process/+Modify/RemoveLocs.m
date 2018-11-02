classdef RemoveLocs<interfaces.DialogProcessor
    methods
        function obj=RemoveLocs(varargin)   
            obj@interfaces.DialogProcessor(varargin{:});  
        end
        
        function out=run(obj,p)
            obj.setPar('undoModule','RemoveLocs');
            notify(obj.P,'backup4undo');
            [~,indroi]=obj.locData.getloc('xnm','position','roi');
            if p.roic.Value==1
            obj.locData.removelocs(indroi);
            else
                obj.locData.removelocs(~indroi);
            end
            obj.locData.regroup;   
            out=[];
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
% pard.textb.object=struct('String','remove locs in ROI','Style','text');
% pard.textb.position=[1,1];
% pard.textb.Width=3;
% pard.textb.object.TooltipString='';

pard.roic.object=struct('String',{{'remove inside ROI','remove outside ROI'}},'Style','popupmenu');
pard.roic.position=[2,1];
pard.roic.Width=2;
pard.roic.object.TooltipString='';
pard.plugininfo.type='ProcessorPlugin';
end