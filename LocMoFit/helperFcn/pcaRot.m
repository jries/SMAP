function out = pcaRot(x,y,varargin)
% Do PCA analysis on the data and determine the rotation angle on xy plane.

p = inputParser;
p.addParameter('xRange_tip',50);
p.addParameter('z',[]);
p.addParameter('output','rotDeg');
p.parse(varargin{:});
xRange_tip = p.Results.xRange_tip;
z = p.Results.z;
output = p.Results.output;

LPC = pca([x y]);
a = LPC(:,1);
b1 = [-1;0];
b2 = [1;0];
ang1 = acos(a'*b1/(norm(a)*norm(b1)));
ang2 = acos(a'*b2/(norm(a)*norm(b2)));
try
    [L1,T] = rotatefactors(LPC);
catch
    warning('Rotation can not be done properly so it is skipped.')
    L1 = [1 0; 0 1];
end
[x1, y1] = rotcoord(x, y, ang1);
[x2, y2] = rotcoord(x, y, ang2);

range1 = range(x1);
range2 = range(x2);

% the minus sign here is for applying the rotatoin to the model.
if range1>range2
    x = x1;
    y = y1;
    rotDeg = rad2deg(ang1);
else
    x = x2;
    y = y2;
    rotDeg = rad2deg(ang2);
end

% detect the position of the most left tip
lInTheRange = x>min(x)&x<min(x)+xRange_tip;
medY = median(y(lInTheRange));
if ~isempty(z)
    medZ = median(z(lInTheRange));
else
    medZ = [];
end
% [mostLeft(1), mostLeft(2)]= rotcoord(mostLeft(1),mostLeft(2), deg2rad(rotDeg));

switch output
    case 'rotDeg'
        out = rotDeg;
    case 'leftTip'
        mostLeft = [min(x) medY medZ];
        out = mostLeft;
    case 'center'
        center = [median(x) median(y) median(z)];
        out = center;
end
end