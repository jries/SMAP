classdef AttachTransformation<interfaces.DialogProcessor

    properties 
    end
    methods
        function obj=AttachTransformation(varargin)   
            obj@interfaces.DialogProcessor(varargin{:}) ;
        end
        
        function out=run(obj,p)
            out=[];
            Tfile=obj.getPar('transformationfile');
            if isempty(Tfile) || contains(Tfile,'temp_T.mat')
                Tfile=strrep(obj.locData.files.file(1).name,'_sml.mat','_T.mat');
            end
            [file,pfad]=uigetfile(Tfile);
            if file
                l=load([pfad file]);
                obj.locData.files.file(1).transformation=l.transformation;

                if p.savesmlfile
                    obj.locData.savelocs(obj.locData.files.file(1).name);
                end
            end
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end




function pard=guidef(obj)
pard.savesmlfile.object=struct('Style','checkbox','String','save SML file','Value',1);
pard.savesmlfile.position=[1,1];
pard.savesmlfile.Width=4;



pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description='';
end




