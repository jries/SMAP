classdef calibrateCMOS<interfaces.DialogProcessor
    methods
        function obj=calibrateCMOS(varargin)      
            obj@interfaces.DialogProcessor(varargin{:}) ;  
        end
        
        function out=run(obj,p)
            out=[];
           
          [file,pfad]=uigetfile('*.tif');
          il=imageloaderAll([pfad  file],[],obj.P);

          mean=double(il.getimage(1));
          M2=0*mean;
          for k=2:il.metadata.numberOfFrames
              image=double(il.getimage(k));
              if isempty(image)
                  break
              end
              delta=image-mean;
              mean=mean+delta/k;
              delta2=image-mean;
              M2=M2+delta.*delta2;
              
          end
          variance=M2/(k-1);
          outputfile=[pfad  strrep(file,'.tif','_var.mat')];
          metadata=il.metadata;
          save(outputfile,'mean','variance','metadata');

             
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function pard=guidef(obj)
% pard.dataselect.object=struct('Style','popupmenu','String','File','Callback',{{@obj.dataselect_callback}});
% pard.dataselect.position=[1,1];
% pard.dataselect.object.TooltipString='sCMOS file localizations';
% pard.dataselect.Width=2;
% pard.tiffselect.object=struct('Style','popupmenu','String','empty');
% pard.tiffselect.position=[2,1];
% pard.tiffselect.object.TooltipString='select Tiff. Use File selector first to populate this menu.';
% pard.tiffselect.Width=2;


pard.plugininfo.type='ProcessorPlugin';
end