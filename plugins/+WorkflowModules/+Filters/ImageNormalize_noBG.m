classdef ImageNormalize_noBG<interfaces.WorkflowModule
%     Performs an Anscombe transform. This converts Poisson noise into
%     Normal distributed noise with unit variance. This can be used to
%     convert the image into a probability map and facilitates segmentation
%     for images with varying background. As ImageNormalize, but for the
%     case when no background is calculated. According to: [1]	U. Koethe,
%     F. Herrmannsdoerfer, I. Kats, and F. A. Hamprecht, SimpleSTORM: a
%     fast, self-calibrating reconstruction algorithm for localization
%     microscopy,HISTOCHEMISTRY AND CELL BIOLOGY, pp. 1-15, Apr. 2014.';
       
    properties
        preview

    end
    methods
        function obj=ImageNormalize_noBG(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=1; 
        end
        function pard=guidef(obj)
            pard.plugininfo.type='WorkflowModule';
            pard.plugininfo.description='As ImageNormalize, but for the case when no background is calculated. Performs an Anscombe transform. This converts Poisson noise into Normal distributed noise with unit variance. This can be used to convert the image into a probability map and facilitates segmentation for images with varying background. According to: [1]	U. Koethe, F. Herrmannsdoerfer, I. Kats, and F. A. Hamprecht, SimpleSTORM: a fast, self-calibrating reconstruction algorithm for localization microscopy,HISTOCHEMISTRY AND CELL BIOLOGY, pp. 1-15, Apr. 2014.';
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
        end
        function prerun(obj,p)          
            obj.preview=obj.getPar('loc_preview');          
        end
        function dato=run(obj,data,p)
            
            if  ~isempty(data.data)
                image=data.data;%get;  
                imnorm=poissonNormalize(image);
                dato=data;%{1}.copy;
                dato.data=imnorm;%set(imnorm);
            else 
                dato=data;
            end     
        end
        

    end
end


function out=poissonNormalize(in)
% out=real(2*sqrt(in+3/8));
in(in<-0.3750)=0;
out=(2*sqrt(in+0.3750));
end
