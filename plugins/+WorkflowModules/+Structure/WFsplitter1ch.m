classdef WFsplitter1ch<WorkflowModules.Structure.WFsplitter
%     This plugin cuts out regions of interest of a defined size around the
%     candidate positions and passes these on to the fitter.
    properties

    end
    methods
        function obj=WFsplitter1ch(varargin)
            obj@WorkflowModules.Structure.WFsplitter(varargin{:})
            obj.inputChannels=1; 
        end
%         function pard=guidef(obj)
%             pard=guidef(obj);
%         end

%         function prerun(obj,p)
% %             splitselection=obj.getSingleGuiParameter('splitWFselection');
%             br=p.splitWFselection.Value;
%             offbranches=1:length(p.splitWFselection.String);
%             offbranches=setdiff(offbranches,br);
%             for b=1:length(offbranches)
%                 br2=offbranches(b);
%                 for k=1:length(obj.modules{br2})
%                     obj.modules{br2}{k}.initialized=true;
%                 end  
%             end 
%          
%         end
%         function nooutput=run(obj,data,p)
%            nooutput=[];
%             switch p.splitWFselection.Value
%                 case 1
%                     obj.output(data,1);
%                 case 2
%                     obj.output(data,2);
%             end
%         end
%         function modelchanged(obj,a,b)
%             splitselection=obj.getSingleGuiParameter('splitWFselection');
%             br=splitselection.Value;
%             offbranches=1:length(splitselection.String);
%             offbranches=setdiff(offbranches,br);
%    
%                 for k=1:length(obj.modules{br})
%                     obj.modules{br}{k}.initialized=false;
%                     guih=obj.modules{br}{k}.guihandles;                    
%                     if isempty(guih)
%                         continue
%                     end
% 
%                     fn=fieldnames(guih);
%                     for l=1:length(fn)
%                         guih.(fn{l}).Visible='on';
%                     end
%                         obj.modules{br}{k}.switchvisibleall;
%                 end
%                 for b=1:length(offbranches)
%                     br2=offbranches(b);
%                     for k=1:length(obj.modules{br2})
%                         guih=obj.modules{br2}{k}.guihandles;
%                         if isempty(guih)
%                             continue
%                         end
%                         fn=fieldnames(guih);
%                         for l=1:length(fn)
%                             guih.(fn{l}).Visible='off';
%                         end
%                     end  
%                 end          
%         end
%         function updateGui(obj)
%            %populate here the splitWFselection menu. 
%             %change visibility of sub-WF
%             for k=1:ceil(length(obj.outputModules))
%                 module=obj.outputModules(k).module;
%                 ind=1;
%                 name{k}='';
%                 while ~contains(module.info.name,'WFcombiner')
%                     obj.modules{k}{ind}=module;
%                     name{k}=[name{k} ',' module.info.name];
%                     ind=ind+1;
%                     module=module.outputModules(1).module;
%                 end
%                 name{k}(1)='';
%             end
%            
%             obj.setGuiParameters(struct('splitWFselection',struct('String',{name},'Value',1))  )  
%             obj.modelchanged;
%         end
    end
end



% function pard=guidef(obj)
% pard.splitWFselection.object=struct('Style','popupmenu','String',{{'1','2'}},'Callback',@obj.modelchanged);
% pard.splitWFselection.position=[1,1];
% pard.splitWFselection.Width=2;
% 
% 
% 
% pard.syncParameters={{'splitWFselection','splitWFselection',{'Value'}}};
% 
% pard.plugininfo.type='WorkflowModule'; 
% pard.plugininfo.description='This plugin cuts out regions of interest of a defined size around the candidate positions and passes these on to the fitter';
% end