classdef saver_settings<interfaces.DialogProcessor
    methods
        function obj=saver_settings(varargin)  
                obj@interfaces.DialogProcessor(varargin{:}) ;
                obj.inputParameters={'filelist_long','mainfile','mainGui','numberOfLayers','sr_layerson'};
        end
        function pard=guidef(obj)
            pard.plugininfo.type='SaverPlugin';
        end
            
        function out=save(obj,p)
            obj.status('save GUI settings')          
            fn=p.filelist_long.selection;          
            [path, file]=fileparts(fn);
            of=[path filesep file '_guisettings.mat'];
            [f,path]=uiputfile(of);
            
            rg=p.mainGui; 
            parameters=rg.saveParameters;
            fileformat.name='guiparameters';
            save([path f],'fileformat','parameters','-v7.3');
            obj.status('save done')
          
        end
        function run(obj,p)
            obj.save(p)
        end
    end
end