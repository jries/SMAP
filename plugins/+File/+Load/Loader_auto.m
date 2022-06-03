classdef Loader_auto<interfaces.DialogProcessor
%     Selects the correct loader based on filename and content
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
            obj.guihandles.updateGuiPar.Value=p.updateGuiPar; %do not overwrite update Guipar settings.
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function out=run(obj,p)
            [f,path]=uigetfile(obj.info.extensions);
            if exist([path f],'file')
                obj.load(p,[path f]);
                initGuiAfterLoad(obj);
                out.file=[f,path];
            else
                out.error='file not found. Cannot be loaded.';
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

pard.restoreROI.object=struct('Style','checkbox','String','restore ROI');
pard.restoreROI.position=[2,1];
pard.restoreROI.Width=2;
pard.restoreROI.TooltipString='Restore ROI if saved with data.';


info.name='auto loader';
info.extensions={'*.mat;*.tif;*.csv;*.hdf5';'*.mat';'*.tif';'*.csv';'*.*'};
info.dialogtitle='select any SMLM position, Tif, csv, settings or workflow file';
pard.plugininfo=info;
pard.plugininfo.type='LoaderPlugin';
pard.plugininfo.description='Selects the correct loader based on filename and content';
end

function loadfile(obj,p,file,mode)            
        [loader,pout]=getloader(obj,mode);
        p=copyfields(p,pout);
        if ~isempty(loader)
        loader.load(p,file);
        obj.notfound=obj.notfound&&loader.notfound;
        end
end

        
function [loader,p]=getloader(obj,mode)  
p=[];
obj.notfound=false;
        switch mode
            case 'tif'
                loadername='Loader_tif';

            case {'sml','fitpos','sites','se'}
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
        switch mode
            case 'hdf5'
                p.importdef.selection='hdf5_Jungmann.txt'; %hack to load Jungmann when automatic loader and hdf5
                p.importdef.Value=8;
        end
end
