function [probability] = nucleationNCoat2(x, y, z, xShift, distanceTop, distanceBottom, ampAll, ampCap, angle, innerRadR, outerRadR, innerRadC, outerRadC, thickness, pxsize, roisize, varargin)
global capmatrix bottommatrix
% Template for a custom PDF
%   Text between double asterisk(**) should be replaced. 
%   See the other .m files in this PDF folder for examples of how to
%   format this type of file 
%   Type the name of your function above, and all your data variables and
%   fitted variables, make the last argument "varargin" if you want to
%   specify a function that can take a tmin and/or tmax

% The metaData variable contains all the information about the PDF
metaData=struct('name','nucleationNCoat',...
            'PDF',  'nucleationNCoat.m',...
            'dataVar','x,y,z',...
            'fitVar', 'xShift,distanceTop,distanceBottom,ampAll,ampCap,angle,innerRadR,outerRadR,innerRadC,outerRadC,thickness,pxsize,roisize',...
            'ub',   'inf,inf,inf,inf,inf,inf,30,70,20,45,60,5,500',...
            'lb',   'inf,inf,inf,inf,inf,inf,30,70,20,45,60,5,500',...
            'guess','inf,inf,inf,inf,inf,inf,30,70,20,45,60,5,500');

% here you can specify the actual functional form if it is too long or 
% complicated to fit into a single line in the metadata above

p.roiSize = 300;

if length(z)==1
    
else
    c = [xShift, distanceTop, distanceBottom, ampAll, ampCap, angle];

    [xrot,yrot] = rotcoord(x,y,-angle*pi/180);

    ytop = round(yrot+distanceTop+p.roiSize/2); ybottom = round(yrot+distanceBottom+p.roiSize/2);
    xall = round(xrot+xShift+p.roiSize/2);
    outOfRange = ytop<1|ytop>p.roiSize|ybottom<1|xall>p.roiSize|xall<1|ytop>p.roiSize;
    
    likTop = zeros(size(ytop));
    likBottom = zeros(size(ybottom));
    
    indtop = sub2ind(size(capmatrix),ytop(~outOfRange), xall(~outOfRange));
    indbottom = sub2ind(size(bottommatrix),ybottom(~outOfRange), xall(~outOfRange));
    
    likTop(~outOfRange) = capmatrix(indtop);
    likBottom(~outOfRange) = bottommatrix(indbottom);
    
    sumLik = sum(likTop)+sum(likBottom);
    
    likTop = likTop/sumLik;
    likBottom = likBottom/sumLik;
    
    probability = (likTop+likBottom++1.1111e-06)/1.1;
end
end