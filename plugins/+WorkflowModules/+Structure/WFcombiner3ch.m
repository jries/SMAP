classdef WFcombiner3ch<interfaces.WorkflowModule
%     This workflow module combines two 
    properties
        loc_ROIsize
        preview
        disppreview
    end
    methods
        function obj=WFcombiner3ch(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=3; 
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
   
            obj.setInputChannels(3,'none');
            %populate here the splitWFselection menu. 
            %change visibility of sub-WF

        end
        function prerun(obj,p)
        end
        function dato=run(obj,data,p)
            dato=data;
%             switch p.splitWFselection.Value
%                 case 1
%                     dato=data{1};
%                 case 2
%                     dato=data{2};
%             end
        end
    end
end



function pard=guidef
% pard.splitWFselection.object=struct('Style','popupmenu','String',{{'1','2'}});
% pard.splitWFselection.position=[1,1];
% pard.splitWFselection.Width=0.8;
% 
% 
% 
% pard.syncParameters={{'splitWFselection','splitWFselection',{'Value'}}};

pard.plugininfo.type='WorkflowModule'; 
pard.plugininfo.description='This plugin cuts out regions of interest of a defined size around the candidate positions and passes these on to the fitter';
end