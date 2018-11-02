classdef DisplayChooser<interfaces.WorkflowModule
    properties
        filterkernel
    end
    methods
        function obj=DisplayChooser(varargin)
            obj@interfaces.WorkflowModule(varargin{:});
            obj.outputParameters={'loc_previewmode'};
            obj.inputChannels=0;
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
%            obj.guihandles.loadmetadata.Callback={@loadmetadata_callback,obj};
%            obj.guihandles.camparbutton.Callback={@camparbutton_callback,obj};
        end
%         function prerun(obj,p)
%       
%         end
%         function run(obj,data,p)
% 
%         end
  
    end
end


function pard=guidef
pard.text.object=struct('Style','text','String','preview mode: ');
pard.text.position=[1,1];
pard.text.Width=1.5;
pard.loc_previewmode.object=struct('Style','popupmenu','String','image-bg|image|norm(image)|bg');
pard.loc_previewmode.position=[2,1];
pard.loc_previewmode.Width=1.5;
pard.loc_previewmode.TooltipString=sprintf('Determine which image to display in Preview mode. Peak finding is performed on norm(image)');
pard.plugininfo.type='WorkflowModule';
pard.plugininfo.description='Select which image to display during preview.';
end