classdef plugin4workflow<interfaces.WorkflowModule
    properties
        subpluginpath
        subplugin
    end
    methods
        function obj=plugin4workflow(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=1;
            obj.isstartmodule=true;
        end

        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
            thisplugin=plugin(obj.subpluginpath{:});
            thisplugin.parent=obj;
            thisplugin.handle=obj.handle;
            thisplugin.attachLocData(obj.locData);
            thisplugin.attachPar(obj.P);
            if isa(thisplugin, 'interfaces.DialogProcessor')
                obj.guiPar.Vrim=0;
                thisplugin.processorgui=false;
            end
            thisplugin.setGuiAppearence(obj.guiPar);
            thisplugin.makeGui;
            
            obj.subplugin=thisplugin;
            obj.children.(thisplugin.pluginpath{end})=thisplugin;

        end
        function prerun(obj,p)
%             p=obj.getGuiParameters;
            module=obj.subplugin;
            if isa(module, 'interfaces.DialogProcessor')
                try
                module.prerun;
                catch
%                     disp('no pre run')
                end
            elseif isa(module, 'interfaces.GuiModuleInterface')
                disp('no dialog processor. implement with run')
            end
        end
        function sethandle(obj,h)
            obj.handle=h;
            obj.subplugin.handle=h;
        end
%         function set.handle(obj,newhandle)
%             obj.handle=newhandle;
%             obj.subplugin.handle=newhandle;
%         end

        function out=run(obj,data,p)
            out=[];
            if ~data.eof
                return
            end
            if isempty(obj.locData.loc)
                disp('no localization data present')
            end
            module=obj.subplugin;
            if isa(module, 'interfaces.DialogProcessor')
                module.processgo;
                try 
                    module.postrun;
                catch
                end
            elseif isa(module, 'interfaces.GuiModuleInterface')
                disp('no dialog processor. implement with run')
            end
            try
                obj.output(data)
            catch err
                
            end
        end

        

    end
end