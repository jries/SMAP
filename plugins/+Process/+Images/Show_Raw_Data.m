classdef Show_Raw_Data<interfaces.DialogProcessor
    methods
        function obj=Show_Raw_Data(varargin)      
            obj@interfaces.DialogProcessor(varargin{:}) ;  
            obj.showresults=true;
        end
        
        function out=run(obj,p)
            out=[];
           
            file=p.dataselect.Value;
          
            rawall=obj.locData.files.file(file).raw;
            for k=length(rawall):-1:1
                imall(:,:,k)=rawall(k).image;
            end
     
            ax=obj.initaxis('image');
            imageslicer(imall,'Parent',ax);
             
            
             
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function pard=guidef(obj)
pard.dataselect.object=struct('Style','popupmenu','String','File');
pard.dataselect.position=[1,1];
pard.dataselect.object.TooltipString='file for target localizations';
pard.dataselect.Width=2;

pard.syncParameters={{'filelist_short','dataselect',{'String'}}};
pard.plugininfo.type='ProcessorPlugin';
end