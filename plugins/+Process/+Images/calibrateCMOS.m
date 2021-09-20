classdef calibrateCMOS<interfaces.DialogProcessor
%     calculates the mean and variancer maps of sCMOS cameras to be used as
%     a camera noise model in the fitter
    methods
        function obj=calibrateCMOS(varargin)      
            obj@interfaces.DialogProcessor(varargin{:}) ;  
        end
        
        function out=run(obj,p)
            out=[];
           
          [file,pfad]=uigetfile('*.tif');
          if ~file
              return
          end
          il=imageloaderAll([pfad  file],[],obj.P);

          t=tic;
          if 0 %naive algorithm
              sumim=double(il.getimage(1));
              sumim2=sumim.^2;
              count=1;
              for k=2:il.metadata.numberOfFrames+1
                  imageh=double(il.getimage(k));
                  if isempty(imageh)
                      break
                  end
                  count=count+1;
                  sumim=sumim+imageh;
                  sumim2=sumim2+imageh.^2;
                  if toc(t)>10
                      disp([num2str(k) ' of ' num2str(il.metadata.numberOfFrames)]);
                      t=tic;
                  end
              end
              mean=sumim/(count);
              variance=sumim2/(count)-mean.^2;
              variance=variance*(count/(count-1));

                  
          else         %Welfords algorithm
              mean=double(il.getimage(1));
              M2=zeros(size(mean),'like',mean); %suggestion Jonas from Birmingham
%               M2=mean.^2;
              count=1;
              for k=2:il.metadata.numberOfFrames+1
                  imageh=double(il.getimage(k));
                  if isempty(imageh)
                      break
                  end
                  count=count+1;
                  delta=imageh-mean; %is this correct?
                  mean=mean+delta/count;
                  delta2=imageh-mean;
                  M2=M2+delta.*delta2;
                  if toc(t)>10
                      disp([num2str(k) ' of ' num2str(il.metadata.numberOfFrames)]);
                      t=tic;
                  end
              end
              variance=M2/(count-1);
          end
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
pard.text.object=struct('Style','text','String','Select tiff stacks containing dark field images of the sCMOS camera.');
pard.text.position=[1,1];
pard.text.Width=4;
% pard.dataselect.object.TooltipString='sCMOS file localizations';
% pard.dataselect.Width=2;
% pard.tiffselect.object=struct('Style','popupmenu','String','empty');
% pard.tiffselect.position=[2,1];
% pard.tiffselect.object.TooltipString='select Tiff. Use File selector first to populate this menu.';
% pard.tiffselect.Width=2;


pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description='calculates the mean and variancer maps of sCMOS cameras to be used as a camera noise model in the fitter';
end