classdef Mathematics<interfaces.DialogProcessor
%     Simple mathematics on localization fields
    methods
        function obj=Mathematics(varargin)      
            obj@interfaces.DialogProcessor(varargin{:}) ;  
        end
        
        function out=run(obj,p) 
            out=[];
            obj.setPar('undoModule','Mathematics');
            notify(obj.P,'backup4undo');
            field=p.assignfield1.selection;
            file=p.dataselect.Value;
            
            locs=obj.locData.getloc({'filenumber','channel',field});
            infile=locs.filenumber==file;
            inchannel=ismember(locs.channel,p.channel);
            inall=infile&inchannel;
            switch p.operator.selection
                case '+'
                 locs.(field)(inall)=locs.(field)(inall)+p.value;
                case '-'
                 locs.(field)(inall)=locs.(field)(inall)-p.value;
                case '*'
                 locs.(field)(inall)=locs.(field)(inall)*p.value;
                case '/'
                 locs.(field)(inall)=locs.(field)(inall)/p.value;
            end
            obj.locData.loc.(field)=locs.(field);
            obj.locData.regroup;
%             
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end

function pard=guidef
pard.texta.object=struct('String','dataset','Style','text');
pard.texta.position=[1,1];

pard.dataselect.object=struct('Style','popupmenu','String','File');
pard.dataselect.position=[2,1];

pard.dataselect.object.TooltipString='choose localization file data set';

pard.textb.object=struct('String','Channel','Style','text');
pard.textb.position=[2,2];
pard.channel.object=struct('Style','edit','String','0 1 2');
pard.channel.position=[2,3];

pard.operator.object=struct('String','+|-|*|/','Style','popupmenu');
pard.operator.position=[3,2];

pard.value.object=struct('Style','edit','String','0');
pard.value.position=[3,3];

pard.assignfield1.object=struct('Style','popupmenu','String','n1|n2');
pard.assignfield1.position=[3,1];

pard.syncParameters={{'filelist_short','dataselect',{'String'}},{'locFields','assignfield1',{'String'}}};

pard.assignfield1.object.TooltipString='choose  field';
pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description='Simple mathematics on localization fields';
end