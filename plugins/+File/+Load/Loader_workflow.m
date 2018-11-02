classdef Loader_workflow<interfaces.DialogProcessor
    methods
        function obj=Loader_workflow(varargin)        
                obj@interfaces.DialogProcessor(varargin{:}) ;
                obj.inputParameters={'mainGui'};
        end
        
        function out=load(obj,p,file,mode)
            if nargin<4
                mode=getfilemode(file);
            end
            loadfile(obj,p,file,mode);
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function run(obj,p)
            [f,p]=uigetfile(obj.info.extensions);
            obj.load(p,[p f]);
            initGuiAfterLoad(obj);
        end
        function clear(file,isadd)
        end
    end
end




function pard=guidef
info.name='workflow loader';
info.extensions={'*.mat';'*.*'};
info.dialogtitle='select workflow file';
pard.plugininfo=info;
pard.plugininfo.type='LoaderPlugin';
end

function loadfile(obj,p,file,mode)            
switch mode
    case 'workflow'
        disp('workflow')
        [~,filename]=fileparts(file);
         module=interfaces.Workflow;
         module.processorgui=false;
         module.handle=figure('MenuBar','none','Toolbar','none','Name',filename);
        module.attachPar(obj.P);
        module.attachLocData(obj.locData);
        p.Vrim=10;
        p.Xrim=10;
        module.setGuiAppearence(p)
        module.makeGui;
        module.guihandles.showresults.Value=1;
        module.load(file);

    otherwise
        disp('file type not recognized')
end
end
        
 