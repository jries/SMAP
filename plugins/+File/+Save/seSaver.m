classdef seSaver<interfaces.DialogProcessor
    methods
        function obj=seSaver(varargin)  
                obj@interfaces.DialogProcessor(varargin{:}) ;
                obj.inputParameters={'filelist_long','mainfile','mainGui','numberOfLayers','sr_layerson','sr_pixrec','layers','sr_image','sr_pos','group_dt','group_dx'};
        end
        
        function out=save(obj,p)
            obj.status('save SE file')
            lastSMLFile = obj.getPar('lastSMLFile');
            defaultFn = replace(lastSMLFile, '_sml', '_se');
            seObj = interfaces.SiteExplorer;
            seObj.sites = obj.locData.SE.sites;
            seObj.cells = obj.locData.SE.cells;
            seObj.files = obj.locData.SE.files;
            uisave('seObj', defaultFn)
            obj.status('save done')
            out = [];
        end
        function pard=guidef(obj)
           pard.plugininfo.type='SaverPlugin';
        end
        function out = run(obj,p)
            obj.save(p)
            out = [];
        end        

    end
end
