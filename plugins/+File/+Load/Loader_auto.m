classdef Loader_auto<interfaces.DialogProcessor
    properties
        notfound=false;
    end
    methods
        function obj=Loader_auto(varargin)        
                obj@interfaces.DialogProcessor(varargin{:}) ;
                obj.inputParameters={'mainGui'};
        end
        
        function out=load(obj,p,file,mode)
            p=copyfields(p,obj.getGuiParameters);
            if nargin<4
                mode=getfilemode(file);
            end
            loadfile(obj,p,file,mode);
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function out=run(obj,p)
            [f,path]=uigetfile(obj.info.extensions);
            obj.load(p,[path f]);
            
            if obj.notfound
                out.error='file not recognized. Cannot be loaded.';
            else
            initGuiAfterLoad(obj);
            end
        end
        function clear(obj,file,isadd)
            mode=getfilemode(file);
            loader=getloader(obj,mode);
            loader.clear(file,isadd);
        end
    end  
end




function pard=guidef

pard.updateGuiPar.object=struct('Style','checkbox','String','load Gui Parameters');
pard.updateGuiPar.position=[1,1];
 pard.updateGuiPar.Width=2;
pard.updateGuiPar.TooltipString='Restore Gui parameters saved with localization data';
info.name='auto loader';
info.extensions={'*.mat;*.tif;*.csv';'*.mat';'*.tif';'*.csv';'*.*'};
info.dialogtitle='select any SMLM position, Tif, csv, settings or workflow file';
pard.plugininfo=info;
pard.plugininfo.type='LoaderPlugin';
end

function loadfile(obj,p,file,mode)            
        loader=getloader(obj,mode);
        if ~isempty(loader)
        loader.load(p,file);
        end
end

        
function loader=getloader(obj,mode)  
obj.notfound=false;
        switch mode
            case 'tif'
                loadername='Loader_tif';

            case {'sml','fitpos','sites'}
                loadername='Loader_sml';
                
            case 'guiparameters'
                loadername='Loader_settings';
       
            case 'workflow'
                loadername='Loader_workflow';
            case {'csv','hdf5','txt'}
                loadername='Loader_csvAndMore';
              
            otherwise
                warning('file type not recognized')
                obj.notfound=true;
                loader=[];
                return
        end
        loader=plugin('File','Load',loadername,[],obj.P);
        loader.attachLocData(obj.locData);
end
