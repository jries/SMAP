classdef ImageNormalize<interfaces.WorkflowModule
%     Performs an Anscombe transform. This converts Poisson noise into
%     Normal distributed noise with unit variance. This can be used to
%     convert the image into a probability map and facilitates segmentation
%     for images with varying background. According to: [1]	U. Koethe, F.
%     Herrmannsdoerfer, I. Kats, and F. A. Hamprecht, SimpleSTORM: a fast,
%     self-calibrating reconstruction algorithm for localization
%     microscopy,HISTOCHEMISTRY AND CELL BIOLOGY, pp. 1-15, Apr. 2014.';
       
    properties
        preview

    end
    methods
        function obj=ImageNormalize(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=2; 
        end
        function pard=guidef(obj)
            pard.plugininfo.type='WorkflowModule';
            pard.plugininfo.description='Performs an Anscombe transform. This converts Poisson noise into Normal distributed noise with unit variance. This can be used to convert the image into a probability map and facilitates segmentation for images with varying background. According to: [1]	U. Koethe, F. Herrmannsdoerfer, I. Kats, and F. A. Hamprecht, SimpleSTORM: a fast, self-calibrating reconstruction algorithm for localization microscopy,HISTOCHEMISTRY AND CELL BIOLOGY, pp. 1-15, Apr. 2014.';
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
            obj.setInputChannels(2,'frame');
        end
        function prerun(obj,p)
            obj.preview=obj.getPar('loc_preview');
        end
        function dato=run(obj,data,p)
            if  ~isempty(data{1}.data)
                image=data{1}.data;%get; 
                bg=data{2}.data;%get;
                indnan=isnan(bg);
                imnorm=poissonNormalize(image)-poissonNormalize(bg);
                imnorm(indnan)=0;
                dato=data{1};%{1}.copy;
                dato.data=imnorm;%set(imnorm);
                if ~obj.preview
%                     obj.output(dato)
                else
                    

                    if data{1}.frame==obj.getPar('loc_previewframe')
                        drawimage(obj,imnorm,image,bg)
                    else
                        dato=[];
                    end
                end
            else 
                dato=data{1};
            end
        end
        

    end
end



function drawimage(obj,imnorm,img,bg)
outputfig=obj.getPar('loc_outputfig');
if ~isvalid(outputfig)
    outputfig=figure(209);
    obj.setPar('loc_outputfig',outputfig);
end

outputfig.Visible='on';
draw=true;
switch obj.getPar('loc_previewmode').Value
    case 1 %image-bg
        imd=img-bg;
    case 2%image
        imd=img;
    case 4 %bg
        imd=bg;
    otherwise 
        draw=false;
end
        
if draw
figure(outputfig)
hold off
imagesc(imd);
colormap jet
colorbar;
axis equal
end
end

function out=poissonNormalize(in)
in(in<-0.3750)=0;
out=(2*sqrt(in+0.3750));
end
