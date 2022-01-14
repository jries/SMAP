function [probability, v1sum, v2sum, nV1, nV2, locxyz,img3d] = NPC3d(x, y, z, xpos, ypos, zpos, distance, rotxz, rotyz, varargin)
    %% only start to modify the codes from nucleationNCoat2_2.m
    p.roiSize = 300;
    axplot = [];
    nargin = length(varargin);
    if nargin >1
        for k = 1:length(varargin(1:2:end))
            ind = k*2-1;
            switch varargin{ind}
                case 'img3d'
                    img3d = varargin{ind+1};
                case 'roisize'
                    roisize = varargin{ind+1};
                case 'SumOffset'
                    SumOffset = varargin{ind+1};
                    sumAllVol = 1+SumOffset;
                    offset = SumOffset/(roisize^2);
                case 'gridX'
                    xcor = varargin{ind+1};
                case 'gridY'
                    ycor = varargin{ind+1};
                case 'gridZ'
                    zcor = varargin{ind+1};
                case 'plot'
                    axplot = varargin{ind+1};
            end
        end
    end


        distanceBottom = -0;
        distanceTop = -(distance-distanceBottom);
        xrot=x-xpos;yrot=y-ypos;zrot=z-zpos;
        [xrot,zrot] = rotcoord(xrot,zrot,-rotxz);
        [yrot,zrot] = rotcoord(yrot,zrot,-rotyz);
%         [xrot,yrot] = rotcoord(xrot,yrot,-rotxy);
        

        ztop = zrot+distanceTop+p.roiSize/2; zbottom = zrot+distanceBottom+p.roiSize/2;       % move the coordinates based on the shift
        xall = xrot+p.roiSize/2;
        yall = yrot+p.roiSize/2;
    if ~isempty(axplot)
        img3d = img3d + shiftImg3DZ(img3d, -distance);
        img3d = shiftImg3DZ(img3d, distance/2);
        imagesc(axplot,[sum(img3d,3),squeeze(sum(img3d,2));squeeze(sum(img3d,1)),0*squeeze(sum(img3d,1))])
        hold(axplot,'on')
        plot(axplot,yall,xall,'k.')
        plot(axplot,(ztop+zbottom)/2+size(img3d,3),xall,'k.')
        
        plot(axplot,(ztop+zbottom)/2,yall+size(img3d,2),'k.')
        hold(axplot,'off')
        title(axplot,['d = ' num2str(distance,3), ', a = ' num2str(rotxz/pi*180,3), ' , b = ' num2str(rotyz/pi*180,3)]);

        return
    else
        outOfRangeTop = yall<1|yall>p.roiSize|xall>p.roiSize|xall<1|ztop>p.roiSize|ztop<1;
        outOfRangeBottom = yall<1|yall>p.roiSize|xall>p.roiSize|xall<1|zbottom>p.roiSize|zbottom<1;

        v1 = zeros(size(zbottom));
        v2 = zeros(size(ztop));

        % interp
        ybottomsub = double(yall(~outOfRangeBottom));
        xbottomsub = double(xall(~outOfRangeBottom));
        zbottomsub = double(zbottom(~outOfRangeBottom));
        ytopsub = double(yall(~outOfRangeTop));
        xtopsub = double(xall(~outOfRangeTop));
        ztopsub = double(ztop(~outOfRangeTop));
        v1(~outOfRangeBottom) = mirt3D_mexinterp(img3d,ybottomsub,xbottomsub,zbottomsub);
%         v1(~outOfRangeBottom) = ba_interp3(img3d,ybottomsub,xbottomsub,zbottomsub,'cubic');
%         interp3(img3d,ybottomsub,xbottomsub,zbottomsub);
%         v2(~outOfRangeTop) = ba_interp3(img3d,ytopsub,xtopsub,ztopsub,'cubic');
        v2(~outOfRangeTop) = mirt3D_mexinterp(img3d,ytopsub,xtopsub,ztopsub);

        nV1 = sum(v1~=0);
        nV2 = sum(v2~=0);
        ampB = 1;
        ampT = nV2/nV1;
        ampTotal = ampB+ampT;
        v1sum = sum(v1); v2sum = sum(v2);
        v1 = v1.*ampB./ampTotal;
        v2 = v2.*ampT./ampTotal;
        v = v1+v2;

        probability = (v + offset)./sumAllVol;
        locxyz = [xall-p.roiSize/2 yall-p.roiSize/2 (ztop+zbottom)/2-p.roiSize/2];
    end
end
function img = shiftImg3DZ(img, shiftZ)
    imgSize = size(img);
    if shiftZ>0
        img = cat(3, zeros([imgSize(1) imgSize(2) ceil(round(shiftZ))]), img(:,:,1:(end-ceil(round(shiftZ)))));
    elseif shiftZ<0
        img = cat(3, img(:,:,(ceil(-round(shiftZ))+1):end), zeros([imgSize(1) imgSize(2) ceil(-round(shiftZ))]));
    else
        return
    end
end