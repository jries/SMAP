classdef imageModel<SMLMModel
    properties (SetObservable, AbortSet)
        pixelSize = 1;
        gaussSigma = 12;
    end
    methods
        function obj = imageModel(img, varargin)
            % parse varargin
            p = inputParser;
            p.addParameter('pixelSize', 2);
            parse(p,varargin{:});
            results = p.Results;
            pixelSize = results.pixelSize;
            
            obj.pixelSize = pixelSize;
            if exist('img')
                obj.img = double(img);
                obj.dimension = length(size(img));
            end
            obj.modelType = 'image';
        end
        function addParent(obj,parent)
            obj.addParent@SMLMModel(parent)
            remVal = rem(parent.roiSize, length(obj.img));
            if remVal==0
                obj.pixelSize = parent.roiSize./length(obj.img);
            else
                warning('Pixel size cannot be detected because the ROI size is not an integer multiple of the image size.')
            end
        end
        function loadListener(obj)
        end
        function img = getImage(obj, mPars, varargin)
            % parse varargin
            p = inputParser;
            p.addParameter('pixelSize', 2);
            p.addParameter('roiSize', 500);
            parse(p,varargin{:});
            results = p.Results;
            pixelSize = results.pixelSize;
            roiSize = results.roiSize;
            
            img = obj.img;
            
            if 1/pixelSize~=1
                if obj.dimension == 2
                    img = imresize(img, 1/(pixelSize/obj.pixelSize));
                else
                    img = imresize3(img, 1/(pixelSize/obj.pixelSize));
                end
            end
        end
        function axs = plot(obj, ringDist)
            fig = figure('Name','imageModel');
            imgSize = size(obj.img);
            if obj.dimension
                [X,Y] = meshgrid(1:imgSize(1),1:imgSize(2)); 
            else
                [X,Y,Z] = meshgrid(1:imgSize(1),1:imgSize(2),1:imgSize(3)); 
            end
            
            F = griddedInterpolant(obj.img,'cubic', 'nearest');
            
            xl=obj.ParentObject.allParsArg.value(1);
            yl=obj.ParentObject.allParsArg.value(2);
            Zl=obj.ParentObject.allParsArg.value(6);
            
            xrot=obj.ParentObject.allParsArg.value(7);
            yrot=obj.ParentObject.allParsArg.value(8);
            
            v1 = F(X(:),Y(:),Z(:));
            newZ = Z-ringDist;
            v2 = F(X(:),Y(:),newZ(:));
            
            
            ind = sub2ind([301 301 301], X(:),Y(:),Z(:));
            vAllImg = zeros([301 301 301]);
            vAllImg(ind) = v1+v2;
            F = griddedInterpolant(vAllImg,'cubic', 'nearest');
            
            
            newX = X-151;
            newY = Y-151;
            newZ = Z-151;
            [newX, newZ] = rotcoord(newX(:)-xl,newZ(:)-Zl,-yrot*pi/180);
            [newY, newZ] = rotcoord(newY(:)-yl,newZ(:),-xrot*pi/180);
            vAll = F(newX+151,newY+151,newZ+151);
            vAllImg = zeros([301 301 301]);
            vAllImg(ind) = vAll;
            
            sideView1 = squeeze(mean(vAllImg,1));
            sideView1 = imadjust(sideView1,prctile(sideView1(:),[50 100]));
            
            sideView2 = squeeze(mean(vAllImg,2));
            sideView2 = imadjust(sideView2,prctile(sideView2(:),[50 100]));
            
            lut = mymakelut('red hot');
            subax1 = axes(fig);
            imagesc(subax1,squeeze(mean(vAllImg,3)'))
            axis equal image
            subplot(2,2,1,subax1);
            subax2 = subplot(2,2,2);
            imagesc(subax2,sideView1)
            axis equal image
            subax3 = subplot(2,2,3);
            imagesc(subax3,sideView2')
            
            % change the colormaps
            colormap(subax1,lut)
            colormap(subax2,lut)
            colormap(subax3,lut)
            
            % rescale the images
            caxis(subax1, [0 0.008]);
            caxis(subax2, [0 1.3]);
            caxis(subax3, [0 1.3]);
            axis equal image
            axs = {subax1, subax2, subax3};
        end
    end
end
