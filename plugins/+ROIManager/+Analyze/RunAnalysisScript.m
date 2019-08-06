classdef RunAnalysisScript<interfaces.DialogProcessor&interfaces.SEProcessor
%     Simple analysis scripts can be saved under +ROIManager/+Anlaysis/scripts. These can be called from the GUI using this plugin.
    properties
        scriptpath='plugins/+ROIManager/+Analyze/scripts/';
    end
    methods
        function obj=RunAnalysisScript(varargin)        
            obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={};
        end
        function initGui(obj)
%             scriptpath='plugins/+ROIManager/+Analyze/scripts/';
            fn=dir([obj.scriptpath '*.m']);
            alln={fn(:).name};
            obj.guihandles.scripts.String=alln;
        end
        
        function out=run(obj,p)  
            [~,funs]=fileparts(p.scripts.selection);
            oldpath=pwd;
            cd(obj.scriptpath)
            funct=str2func(funs);
            try
            funct()
            catch err
                warning err
            end
            cd(oldpath)

            out=[];
           
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end


function pard=guidef(obj)

pard.scripts.object=struct('String',{{''}},'Style','listbox');
pard.scripts.position=[5,1];
pard.scripts.Width=4;
pard.scripts.Height=5;
pard.plugininfo.description='Simple analysis scripts can be saved under +ROIManager/+Anlaysis/scripts. These can be called from the GUI using this plugin.';
pard.plugininfo.type='ROI_Analyze';
end