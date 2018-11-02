function tform = findAffineTransformZ(Xr,zt)
%
% For an affine transformation:
%
%
%                     [ A D 0 ]
% [u v 1] = [x y 1] * [ B E 0 ]
%                     [ C F 1 ]
%
% There are 6 unknowns: A,B,C,D,E,F
%
% Another way to write this is:
%
%                   [ A D ]
% [u v] = [x y 1] * [ B E ]
%                   [ C F ]
%
% Rewriting the above matrix equation:
% U = X * T, where T = reshape([A B C D E F],3,2)
%
% With 3 or more correspondence points we can solve for T,
% T = X\U which gives us the first 2 columns of T, and
% we know the third column must be [0 0 1]'.

% [uv,normMatrix1] = images.geotrans.internal.normalizeControlPoints(uv);
% [xy,normMatrix2] = images.geotrans.internal.normalizeControlPoints(xy);

minRequiredNonCollinearPairs = 3;
M = size(Xr,1);
X = [Xr ones(M,1)];

% just solve for the first two columns of T
U = zt;

% We know that X * T = U
if rank(X)>=minRequiredNonCollinearPairs

    Tinv = X \ U;
else
    error(message('images:geotrans:requiredNonCollinearPoints', minRequiredNonCollinearPairs, 'affine'))    
end

% add third column
Tinvb=eye(4);
Tinvb(:,3)=Tinv;
% Tinv(:,3) = [0 0 1]';

% Tinv = normMatrix2 \ (Tinv * normMatrix1);

T = inv(Tinvb);
% T(:,3) = [0 0 1]';

tform = affine3d(Tinvb);