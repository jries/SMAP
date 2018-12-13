classdef fibrilFinder<interfaces.DialogProcessor&interfaces.SEProcessor
    % This is a plugin in development. Public has no access to the 
    % run_ functions called in this plugin. For internal users,
    % 'fibrilKymograph' is required.
    methods
        function obj=fibrilFinder(varargin) 
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={};
        end
        function out=run(obj,p)
            out=runintern(obj,p);
         
        end
        
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function pard = guidef(obj)
    pard.t1.object = struct('Style','text','String','test');
    pard.t1.position = [1,1];
    pard.t1.Width = 1;
    
    pard.plugininfo.type='ROI_Segment';
end

function out=runintern(obj,p)
    out = [];
    p.pxSize = 7;
    p.minArea = 5000;
    files = obj.locData.SE.files;
    for n = 1:obj.locData.SE.numberOfFiles
        fileID = files(n).ID;
        allLocs = obj.locData.getloc({'xnm', 'ynm', 'frame'}, 'grouping', 'grouped', 'layer',1, 'filenumber',fileID);
        for cn = 1:obj.locData.SE.numberOfCells
            % for each file, get its locs
            locs = allLocs;
            indInCell = obj.locData.SE.cells(cn).pos(1) + obj.locData.getPar('se_cellfov')/2 & locs.xnm >= obj.locData.SE.cells(cn).pos(1) - obj.locData.getPar('se_cellfov')/2;
            indInCell = indInCell & obj.locData.SE.cells(cn).pos(2) + obj.locData.getPar('se_cellfov')/2 & locs.xnm >= obj.locData.SE.cells(cn).pos(2) - obj.locData.getPar('se_cellfov')/2;
            locs.xnm = allLocs.xnm(indInCell,1);
            locs.ynm = allLocs.ynm(indInCell,1);
            locs.frame = allLocs.frame(indInCell,1);
            
            % rendering the image
            minx = min(locs.xnm);
            miny = min(locs.ynm);
            xrange = p.pxSize*floor(min(locs.xnm)/p.pxSize):p.pxSize:p.pxSize*ceil(max(locs.xnm)/p.pxSize);
            yrange = p.pxSize*floor(min(locs.ynm)/p.pxSize):p.pxSize:p.pxSize*ceil(max(locs.ynm)/p.pxSize);
            img = histcounts2(locs.ynm,locs.xnm, yrange,xrange);
            
            % filtering the image and removing noise
            h = fspecial('gaussian',round(p.pxSize*1.2),round(p.pxSize*1.2)/2);
            img2 = filter2(h,img);
            SE = strel('disk',10);
%             img2 = imtophat(img2,SE);
%             img2 = medfilt2(img2);
            
            % filtering the images by the pattern of stripes
                % creat patterns
            ele1 = zeros([2 15]);
            ele2 = ones([2 15]);
            pattern = repmat([ele1;ele2],[5 1]);
            pattern = pattern/sum(pattern(:));
            pattern2 = imrotate(pattern, 45, 'crop');
            pattern3 = imrotate(pattern, -45, 'crop');
                % filtering
            imgFilt = filter2(pattern2, img2);
            imgFilt = filter2(pattern3, imgFilt);
            figure; imagesc(imgFilt);
            
            
            % remove noise
            imgFilt = medfilt2(imgFilt);
            imgFilt = imtophat(imgFilt,strel('disk',10));
%             imgbw = imerode(imgbw, strel('disk', 5));
%             imgFilt = medfilt2(imgFilt);

            % create mask
            imgbw = imgFilt > 0.01;
%             imgbw = imopen(imgbw,2);
            imgbw = bwareafilt(imgbw,[1000/p.pxSize^2*100 inf]);
            figure; imagesc(imgbw)
            
            % detect ridges
            Df = bwdist(~imgbw);
            h = fspecial('gaussian', 7,4);
            DfFilt = filter2(h, Df);
            [gx, gy] = gradient(double(-DfFilt));
            [gxx, gxy] = gradient(gx);
            [gxy, gyy] = gradient(gy);
            
            eigVal1 = (gxx + gyy + sqrt((gxx - gyy).^2 + 4*gxy.^2))/2;
%             eigVal2 = (gxx + gyy - sqrt((gxx - gyy).^2 + 4*gxy.^2))/2;
%             Get Prominence and local maxima of distance to nearst
%             non-zero pixels
            ridge = eigVal1 > 0.2;
            ridge = double(ridge);
            imgbw = imopen(imgbw, strel('disk',5));
            figure; imagesc(imgbw + (imgbw&ridge))
            axis('equal')
            % the standard skeletonization:
%             skelimg = bwmorph(ridge,'thin',inf);
% 
%             mn = bwmorph(skelimg,'branchpoints');
%             mn = imclose(mn, strel('disk',3));
%             [row, column] = find(mn);
%             branchpts = [row column];
% 
% 
%             Endimg = bwmorph(skelimg,'endpoints');
%             Endimg = imclose(Endimg, strel('disk',20));
%             Endimg = bwmorph(Endimg,'shrink',inf);
%             [row,column] = find(Endimg);
%             Endpts = [row column];
%             
% 
%             n = size(Endpts,1);
%             Cntrpts = zeros(n,2);
%             paths = zeros(size(imgbw));
%             for ii = 1:n
%                 % compute end & branch points geodesic distance transform
%                 dEnd = bwdistgeodesic(imgbw, Endpts(ii,2), Endpts(ii,1), 'quasi-euclidean');
%                 [~,closestBranchIdx] = min(dEnd(mn));
%                 dStart = bwdistgeodesic(imgbw, branchpts(closestBranchIdx,2), branchpts(closestBranchIdx,1), 'quasi-euclidean');
%                 D = dStart + dEnd;
%                 D = round(D * 8) / 8;
%                 D(isnan(D)) = inf;
%                 figure(954485); imagesc(imextendedmin(D,5))
%                 pause(0.001)
%             end
            
            % detect fibril
            imgbwfilt = bwareafilt(logical(imgbw), [p.minArea/p.pxSize inf]);
            fibril = bwconncomp(imgbwfilt);
            
            preFibrilCor = zeros([fibril.NumObjects 2]);
            xcor = cell([fibril.NumObjects 1]);
            ycor = cell([fibril.NumObjects 1]);
            for k = 1:fibril.NumObjects
                [preFibrilCor(k,:), xcor{k}, ycor{k}]= getMidPoints(fibril.ImageSize,fibril.PixelIdxList{k});
            end
            fibrilCor = preFibrilCor*p.pxSize;

            indRidge = find(ridge(:));                  % determine the coordinates of ridges
            
            % transform locs to pixel positions
            oneX = floor((locs.ynm-miny)/p.pxSize);
            oneX(oneX==0)=1;                            % to solve the zero generated by floor()
            oneY = floor((locs.xnm-minx)/p.pxSize);
            oneY(oneY==0)=1;                            % to solve the zero generated by floor()
%                 get the indice of locs in the image
            inImg = sub2ind(size(imgbw), oneX, oneY);
            
            % get site-specific info.
            for k = 1:fibril.NumObjects
                segmentation = [];
                indOneFibril = ismember(inImg,fibril.PixelIdxList{k});
                indOneRidge = ismember(indRidge,fibril.PixelIdxList{k});
                [subOneRidgey, subOneRidgex] = ind2sub(size(imgbw), indRidge(indOneRidge));
                segmentation.indOneFibril = indOneFibril;
                segmentation.subOneRidgey = subOneRidgey;
                segmentation.subOneRidgex = subOneRidgex;
                segmentation.minx = minx;
                segmentation.miny = miny;
                segmentation.pxSize= p.pxSize;
                currentsite=interfaces.SEsites;
                currentsite.pos = [fibrilCor(k,2)+min(locs.xnm) fibrilCor(k,1)+min(locs.ynm)];
                currentsite.evaluation.segmentation = segmentation;
                currentsite.ID=0;
                currentsite.info.filenumber=fileID;
                currentsite.info.cell=cn;
                obj.locData.SE.addSite(currentsite);
            end   
        end
    end
end

function [sub, a, b] = getMidPoints(imageSize, idxlist)
    [a,b] = ind2sub(imageSize, idxlist);
    px = [a b];
    sub = (min(px)+max(px))/2;
end