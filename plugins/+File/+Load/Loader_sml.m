classdef Loader_sml<interfaces.DialogProcessor
%     Loads  localization data in the SMAP proprietary format '_sml.mat'
    methods
        function obj=Loader_sml(varargin)        
                obj@interfaces.DialogProcessor(varargin{:}) ;
                obj.inputParameters={'mainGui'};
        end
        
        function out=load(obj,p,file,mode)
            if nargin<4
                mode=getfilemode(file);
            end
            if isempty(p)
                p.updateGuiPar=false;
            end
            loadfile(obj,p,file,mode);
        end
        function run(obj,p)
            [f,p]=uigetfile(obj.info.extensions);
            obj.load(p,[p f]);
            initGuiAfterLoad(obj);
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function clear(obj,file,isadd)
            if nargin<3
                isadd=false;
            end
            if isadd 
                obj.locData.clear('filter');
            else
                obj.locData.clear;
            end
        end        
    end
end




function pard=guidef
info.name='SML loader';
info.extensions={'*.mat;*.*'};
info.dialogtitle='select any SMLM file (_sml or _fitpos)';
pard.plugininfo=info;
pard.plugininfo.type='LoaderPlugin';
pard.plugininfo.description='Loads  localization data in the SMAP proprietary format _sml.mat';
end

