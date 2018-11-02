classdef autoLocSaver<interfaces.WorkflowModule;
    properties
        
    end
    methods
       function obj=autoLocSaver(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=1; 
            obj.isstartmodule=true;
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
        end
        function prerun(obj,p)
            
        end
        function output=run(obj,data,p)
            output=[];
            fn=obj.getPar('lastSMLFile');
            if isempty(fn)
                fn=obj.locData.files.file(1).name;
            end
            fn=strrep(fn,'_sml',['_' p.postfix '_sml']);
            try
                obj.locData.save(findunusedfilename(fn));
            catch err
                err
            end
           
        end
   
    end
end



function pard=guidef(obj)
pard.t1.object=struct('Style','text','String','add to filename:');
pard.t1.position=[5,1];
pard.t1.Width=4;

pard.postfix.object=struct('Style','edit','String','wf');
pard.postfix.position=[6,1];
pard.postfix.Width=0.5;

end