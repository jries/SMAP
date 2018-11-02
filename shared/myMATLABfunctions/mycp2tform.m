function [trans,uv,xy,uv_dev,xy_dev] = cp2tform(varargin)
%CP2TFORM Infer spatial transformation from control point pairs.
%   CP2TFORM takes pairs of control points and uses them to infer a
%   spatial transformation. 
%
%   TFORM = CP2TFORM(INPUT_POINTS,BASE_POINTS,TRANSFORMTYPE) returns a TFORM
%   structure containing a spatial transformation. INPUT_POINTS is an M-by-2
%   double matrix containing the X and Y coordinates of control points in
%   the image you want to transform. BASE_POINTS is an M-by-2 double matrix
%   containing the X and Y coordinates of control points in the base
%   image. TRANSFORMTYPE can be 'nonreflective similarity', 'similarity',
%   'affine', 'projective', 'polynomial', 'piecewise linear' or 'lwm'. See the
%   reference page for CP2TFORM for information about choosing TRANSFORMTYPE.
%
%   TFORM = CP2TFORM(CPSTRUCT,TRANSFORMTYPE) works on a CPSTRUCT structure
%   that contains the control point matrices for the input and base
%   images. The Control Point Selection Tool, CPSELECT, creates the
%   CPSTRUCT.
%
%   [TFORM,INPUT_POINTS,BASE_POINTS] = CP2TFORM(CPSTRUCT,...) returns the
%   control points that were actually used in INPUT_POINTS, and
%   BASE_POINTS. Unmatched and predicted points are not used. See
%   CPSTRUCT2PAIRS.
%
%   TFORM = CP2TFORM(INPUT_POINTS,BASE_POINTS,'polynomial',ORDER) 
%   ORDER specifies the order of polynomials to use. ORDER can be 2, 3, 4.
%   If you omit ORDER, it defaults to 3.
%
%   TFORM = CP2TFORM(CPSTRUCT,'polynomial',ORDER) works on a CPSTRUCT
%   structure.
%
%   TFORM = CP2TFORM(INPUT_POINTS,BASE_POINTS,'piecewise linear') Creates a
%   delaunay triangulation of the base control points, and maps
%   corresponding input control points to the base control points. The
%   mapping is linear (affine) for each triangle, and continuous across the
%   control points, but not continuously differentiable as each triangle has
%   its own mapping.
%
%   TFORM = CP2TFORM(CPSTRUCT,'piecewise linear') works on a CPSTRUCT
%   structure.
%
%   TFORM = CP2TFORM(INPUT_POINTS,BASE_POINTS,'lwm',N) The local weighted
%   mean (lwm) method creates a mapping, by inferring a polynomial at each
%   control point using neighboring control points. The mapping at any
%   location depends on a weighted average of these polynomials.  You can
%   optionally specify the number of points, N, used to infer each
%   polynomial. The N closest points are used to infer a polynomial of order
%   2 for each control point pair. If you omit N, it defaults to 12. N can
%   be as small as 6, BUT making N small risks generating ill-conditioned
%   polynomials.
%
%   TFORM = CP2TFORM(CPSTRUCT,'lwm',N) works on a CPSTRUCT structure.
%
%   [TFORM,INPUT_POINTS,BASE_POINTS,INPUT_POINTS_BAD,BASE_POINTS_BAD] = ...
%        CP2TFORM(INPUT_POINTS,BASE_POINTS,'piecewise linear') 
%   returns the control points that were actually used in INPUT_POINTS and
%   BASE_POINTS, and the control points that were eliminated because they
%   were middle vertices of degenerate fold-over triangles in
%   INPUT_POINTS_BAD and BASE_POINTS_BAD.
%
%   [TFORM,INPUT_POINTS,BASE_POINTS,INPUT_POINTS_BAD,BASE_POINTS_BAD] = ...
%        CP2TFORM(CPSTRUCT,'piecewise linear') works on a CPSTRUCT structure.
%
%   TRANSFORMTYPE
%   -------------
%   CP2TFORM requires a minimum number of control point pairs to infer a
%   TFORM structure of each TRANSFORMTYPE:
%
%       TRANSFORMTYPE         MINIMUM NUMBER OF PAIRS
%       -------------         -----------------------
%       'nonreflective similarity'       2 
%       'similarity'                     3 
%       'affine'                         3 
%       'projective'                     4 
%       'polynomial' (ORDER=2)           6
%       'polynomial' (ORDER=3)          10
%       'polynomial' (ORDER=4)          15
%       'piecewise linear'               4
%       'lwm'                            6
%      
%   When TRANSFORMTYPE is 'nonreflective similarity', 'similarity', 'affine',
%   'projective', or 'polynomial', and INPUT_POINTS and BASE_POINTS (or
%   CPSTRUCT) have the minimum number of control points needed for a particular
%   transformation, the coefficients are found exactly. If INPUT_POINTS and
%   BASE_POINTS have more than the minimum, a least squares solution is
%   found. See MLDIVIDE.
%
%   Note
%   ----
%   When either INPUT_POINTS or BASE_POINTS has a large offset with
%   respect to their origin (relative to range of values that it spans), the
%   points are shifted to center their bounding box on the origin before
%   fitting a TFORM structure.  This enhances numerical stability and
%   is handled transparently by wrapping the origin-centered TFORM within a
%   custom TFORM that automatically applies and undoes the coordinate shift
%   as needed. This means that fields(T) may give different results for
%   different coordinate inputs, even for the same TRANSFORMTYPE.
%
%   Example
%   -------
%   I = checkerboard;
%   J = imrotate(I,30);
%   base_points = [11 11; 41 71];
%   input_points = [14 44; 70 81];
%   cpselect(J,I,input_points,base_points);
%
%   t = cp2tform(input_points,base_points,'nonreflective similarity');
%
%   % Recover angle and scale by checking how a unit vector 
%   % parallel to the x-axis is rotated and stretched. 
%   u = [0 1]; 
%   v = [0 0]; 
%   [x, y] = tformfwd(t, u, v); 
%   dx = x(2) - x(1); 
%   dy = y(2) - y(1); 
%   angle = (180/pi) * atan2(dy, dx) 
%   scale = 1 / sqrt(dx^2 + dy^2) 
%
%  See also CPSELECT, CPCORR, CPSTRUCT2PAIRS, IMTRANSFORM, TFORMFWD, TFORMINV.

%   Copyright 1993-2009 The MathWorks, Inc. 
%   $Revision: 1.10.4.17 $  $Date: 2010/11/17 11:23:34 $

% Note: 'linear conformal' was deprecated in R2008a in favor of 
% 'nonreflective similarity'. Both will still work and give the same 
% result as each other.

[uv, xy, method, options] = ParseInputs(varargin{:});

% initialize deviation matrices
xy_dev = [];
uv_dev = [];

% Assign function according to method and 
% set K = number of control point pairs needed. 
switch method
  case 'nonreflective similarity'
    findT_fcn = @findNonreflectiveSimilarity;
  case 'similarity'
    findT_fcn = @findSimilarity;
  case 'affine'
    findT_fcn = @findAffineTransform;
  case 'projective'
    findT_fcn = @findProjectiveTransform;
  case 'polynomial'
    findT_fcn = @findPolynomialTransform;
  case 'piecewise linear'
    findT_fcn = @findPiecewiseLinear;
  case 'lwm'
    findT_fcn = @findLWM;
  otherwise
    eid = sprintf('Images:%s:internalProblem', mfilename);
    error(eid,...
          'Function %s has an internal problem: unrecognized method.',...
          mfilename);
end    

% error if user enters too few control point pairs
M = size(uv,1);
if M<options.K
    [msg, eid] = CountError(options.K,method);
    error(eid,msg);
end

% get offsets to apply to before/after spatial transformation
uvShift = getShift(uv);
xyShift = getShift(xy);
needToShift = any([uvShift xyShift] ~= 0);

if ~needToShift
    % infer transform
    [trans, output] = findT_fcn(uv,xy,options);
else
    % infer transform for shifted data
    [tshifted, output] = findT_fcn(applyShift(uv,uvShift),...
                                   applyShift(xy,xyShift),options);

    % construct custom tform with tshifted between forward and inverse shifts
    tdata = struct('uvShift',uvShift,'xyShift',xyShift,'tshifted',tshifted);
    trans = maketform('custom',2,2,@fwd,@inverse,tdata);
end

if strcmp(method,'piecewise linear')
    uv = undoShift(output.uv,uvShift);  % No-ops if needToShift
    xy = undoShift(output.xy,xyShift);  % is false.
    uv_dev = output.uv_dev;
    xy_dev = output.xy_dev;
end

%-------------------------------
function shift = getShift(points)
tol = 1e+3;
minPoints = min(points);
maxPoints = max(points);
center = (minPoints + maxPoints) / 2;
span = maxPoints - minPoints;
if (span(1) > 0 && abs(center(1))/span(1) > tol) ||...
   (span(2) > 0 && abs(center(2))/span(2) > tol)
    shift = center;
else
    shift = [0 0];
end

%-------------------------------
function shiftedPoints = applyShift(points,shift)
shiftedPoints = bsxfun(@minus, points, shift);

%-------------------------------
function points = undoShift(shiftedPoints,shift)
points = bsxfun(@plus, shiftedPoints, shift);

%-------------------------------
function x = fwd(u,t)
x = undoShift(tformfwd(applyShift(u,t.tdata.uvShift),t.tdata.tshifted),...
                 t.tdata.xyShift);

%-------------------------------
function u = inverse(x,t)
u = undoShift(tforminv(applyShift(x,t.tdata.xyShift),t.tdata.tshifted),...
                 t.tdata.uvShift);

%-------------------------------
% Function  findAffineTransform
%
function [trans, output] = findAffineTransform(uv,xy,options)
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

K = options.K;
M = size(xy,1);
X = [xy ones(M,1)];

% just solve for the first two columns of T
U = uv;

% We know that X * T = U
if rank(X)>=K
    Tinv = X \ U;
else
    [msg, eid] = RankError(K,'affine');
    error(eid,msg);
end

% add third column
Tinv(:,3) = [0 0 1]';

T = inv(Tinv);
T(:,3) = [0 0 1]';

trans = maketform('affine', T);
output = [];

%-------------------------------
% Function  findProjectiveTransform
%
function [trans, output] = findProjectiveTransform(uv,xy,options)
%
% For a projective transformation:
%
% u = (Ax + By + C)/(Gx + Hy + I)
% v = (Dx + Ey + F)/(Gx + Hy + I)
%
% Assume I = 1, multiply both equations, by denominator:
%
% u = [x y 1 0 0 0 -ux -uy] * [A B C D E F G H]'
% v = [0 0 0 x y 1 -vx -vy] * [A B C D E F G H]'
%
% With 4 or more correspondence points we can combine the u equations and
% the v equations for one linear system to solve for [A B C D E F G H]:
%
% [ u1  ] = [ x1  y1  1  0   0   0  -u1*x1  -u1*y1 ] * [A]
% [ u2  ] = [ x2  y2  1  0   0   0  -u2*x2  -u2*y2 ]   [B]
% [ u3  ] = [ x3  y3  1  0   0   0  -u3*x3  -u3*y3 ]   [C]
% [ u1  ] = [ x4  y4  1  0   0   0  -u4*x4  -u4*y4 ]   [D]
% [ ... ]   [ ...                                  ]   [E]
% [ un  ] = [ xn  yn  1  0   0   0  -un*xn  -un*yn ]   [F]
% [ v1  ] = [ 0   0   0  x1  y1  1  -v1*x1  -v1*y1 ]   [G]
% [ v2  ] = [ 0   0   0  x2  y2  1  -v2*x2  -v2*y2 ]   [H]
% [ v3  ] = [ 0   0   0  x3  y3  1  -v3*x3  -v3*y3 ]
% [ v4  ] = [ 0   0   0  x4  y4  1  -v4*x4  -v4*y4 ]
% [ ... ]   [ ...                                  ]  
% [ vn  ] = [ 0   0   0  xn  yn  1  -vn*xn  -vn*yn ]
%
% Or rewriting the above matrix equation:
% U = X * Tvec, where Tvec = [A B C D E F G H]'
% so Tvec = X\U.
%

K = options.K;
M = size(xy,1);
x = xy(:,1);
y = xy(:,2);
vec_1 = ones(M,1);
vec_0 = zeros(M,1);
u = uv(:,1);
v = uv(:,2);

U = [u; v];

X = [x      y      vec_1  vec_0  vec_0  vec_0  -u.*x  -u.*y;
     vec_0  vec_0  vec_0  x      y      vec_1  -v.*x  -v.*y  ];

% We know that X * Tvec = U
if rank(X) >= 2*K 
    Tvec = X \ U;    
else
    [msg, eid] = RankError(K,'projective');
    error(eid,msg);
end

% We assumed I = 1;
Tvec(9) = 1;

Tinv = reshape(Tvec,3,3);
T = inv(Tinv);

trans = maketform('projective', T);
output = [];

%--------------------------------------
% Function  findNonreflectiveSimilarity
%
function [trans, output] = findNonreflectiveSimilarity(uv,xy,options)
%
% For a nonreflective similarity:
%
% let sc = s*cos(theta)
% let ss = s*sin(theta)
%
%                   [ sc -ss
% [u v] = [x y 1] *   ss  sc
%                     tx  ty]
%
% There are 4 unknowns: sc,ss,tx,ty.
%
% Another way to write this is:
%
% u = [x y 1 0] * [sc
%                  ss
%                  tx
%                  ty]
%
% v = [y -x 0 1] * [sc
%                   ss
%                   tx
%                   ty]
%
% With 2 or more correspondence points we can combine the u equations and
% the v equations for one linear system to solve for sc,ss,tx,ty.
%
% [ u1  ] = [ x1  y1  1  0 ] * [sc]
% [ u2  ]   [ x2  y2  1  0 ]   [ss]
% [ ... ]   [ ...          ]   [tx]
% [ un  ]   [ xn  yn  1  0 ]   [ty]
% [ v1  ]   [ y1 -x1  0  1 ]   
% [ v2  ]   [ y2 -x2  0  1 ]    
% [ ... ]   [ ...          ]
% [ vn  ]   [ yn -xn  0  1 ]
%
% Or rewriting the above matrix equation:
% U = X * r, where r = [sc ss tx ty]'
% so r = X\U.
%

K = options.K;
M = size(xy,1);
x = xy(:,1);
y = xy(:,2);
X = [x   y  ones(M,1)   zeros(M,1);
     y  -x  zeros(M,1)  ones(M,1)  ];

u = uv(:,1);
v = uv(:,2);
U = [u; v];

% We know that X * r = U
if rank(X) >= 2*K 
    r = X \ U;    
else
    msg = 'At least 2 unique points needed to infer nonreflective similarity transform.';
    eid = sprintf('Images:%s:twoUniquePointsReq',mfilename);
    error(eid,msg);
end

sc = r(1);
ss = r(2);
tx = r(3);
ty = r(4);

Tinv = [sc -ss 0;
        ss  sc 0;
        tx  ty 1];

T = inv(Tinv);
T(:,3) = [0 0 1]';

trans = maketform('affine', T);
output = [];

%-------------------------
% Function  findSimilarity
%
function [trans, output] = findSimilarity(uv,xy,options)
%
% The similarities are a superset of the nonreflective similarities as they may
% also include reflection.
%
% let sc = s*cos(theta)
% let ss = s*sin(theta)
%
%                   [ sc -ss
% [u v] = [x y 1] *   ss  sc
%                     tx  ty]
%
%          OR
%
%                   [ sc  ss
% [u v] = [x y 1] *   ss -sc
%                     tx  ty]
%
% Algorithm:
% 1) Solve for trans1, a nonreflective similarity.
% 2) Reflect the xy data across the Y-axis, 
%    and solve for trans2r, also a nonreflective similarity.
% 3) Transform trans2r to trans2, undoing the reflection done in step 2.
% 4) Use TFORMFWD to transform uv using both trans1 and trans2, 
%    and compare the results, returning the transformation corresponding 
%    to the smaller L2 norm.

% Need to reset options.K to prepare for calls to findNonreflectiveSimilarity.
% This is safe because we already checked that there are enough point pairs.
options.K = 2;

% Solve for trans1
[trans1, output] = findNonreflectiveSimilarity(uv,xy,options);


% Solve for trans2

% manually reflect the xy data across the Y-axis
xyR = xy;
xyR(:,1) = -1*xyR(:,1);

trans2r  = findNonreflectiveSimilarity(uv,xyR,options);

% manually reflect the tform to undo the reflection done on xyR
TreflectY = [-1  0  0;
              0  1  0;
              0  0  1];
trans2 = maketform('affine', trans2r.tdata.T * TreflectY);


% Figure out if trans1 or trans2 is better
xy1 = tformfwd(trans1,uv);
norm1 = norm(xy1-xy);

xy2 = tformfwd(trans2,uv);
norm2 = norm(xy2-xy);

if norm1 <= norm2
    trans = trans1;
else
    trans = trans2;
end

%-------------------------------
% Function  findPolynomialTransform
%
function [trans, output] = findPolynomialTransform(uv,xy,options)
%
% For a polynomial transformation: 
%
% u = X*A, v = X*B, solve for A and B:
%     A = X\u;
%     B = X\v;   
% 
% The matrix X depends on the order of the polynomial.
% X will be M-by-K, where K = (order+1)*(order+2)/2;
% so A and B will be vectors of length K.
%
%   order = 2   
%     X = [ones(M,1),  x,  y,  x.*y,  x.^2,  y.^2];
%     so X is an M-by-6 matrix 
%
%   order = 3
%     X = [ones(M,1),  x,  y,  x.*y,  x.^2,  y.^2, ...
%          (x.^2).*y,  (y.^2).*x,  x.^3,  y.^3];
%     so X is an M-by-10 matrix 
%
%   order = 4
%     X = [ones(M,1),  x,  y,  x.*y,  x.^2,  y.^2, ...
%          (x.^2).*y,  (y.^2).*x,  x.^3,  y.^3, ...
%          (x.^3).*y,  (x.^2).*(y.^2),  x.*(y.^3),  x.^4,  y.^4];
%     so X is an M-by-15 matrix 
%
%

K = options.K;
order = options.order;

u = uv(:,1);
v = uv(:,2);

X = getTerms(order,xy);

if rank(X)>=K
    % u = X*A, v = X*B, solve for A and B:
    A = X\u;
    B = X\v;   
else
    [msg, eid] = RankError(K,'polynomial');
    error(eid,msg);
end

trans = maketform('custom',2,2,[],@inv_polynomial,[A B]);
output = [];

%-------------------------------
%
%
function uv = inv_polynomial(xy,t)

% xy must be M-by-2
% t.tdata must be nterms-by-2, where nterms is 6, 10, or 15.

if size(xy,2) ~= 2
    msg = 'XY must have two columns.';
    eid = sprintf('Images:%s:xyMustHave2Cols',mfilename);
    error(eid,msg);
end

if size(xy,1) < 1
    msg = 'XY must have at least one row.';
    eid = sprintf('Images:%s:xyMustHaveOneRow',mfilename);
    error(eid,msg);
end

if size(t.tdata,2) ~= 2
    msg = 'TDATA must have two columns.';
    eid = sprintf('Images:%s:tdataMustHave2Cols',mfilename);
    error(eid,msg);
end

nterms = size(t.tdata,1);

switch nterms
  case 6   
    order = 2;
    
  case 10
    order = 3;
    
  case 15
    order = 4;
    
  otherwise
    msg = 'TDATA must have 6, 10, or 15 rows.';
    eid = sprintf('Images:%s:tdataInvalidNumOfRows',mfilename);
    error(eid,msg);
end

X = getTerms(order,xy);

uv = X * t.tdata;

%-------------------------------
%
%
function X = getTerms(order,xy)

M = size(xy,1);
x = xy(:,1);
y = xy(:,2);

switch order
  case 2   
    X = [ones(M,1),  x,  y,  x.*y,  x.^2,  y.^2];
    
  case 3
    X = [ones(M,1),  x,  y,  x.*y,  x.^2,  y.^2, ...
         (x.^2).*y,  (y.^2).*x,  x.^3,  y.^3];
    
  case 4
    X = [ones(M,1),  x,  y,  x.*y,  x.^2,  y.^2, ...
         (x.^2).*y,  (y.^2).*x,  x.^3,  y.^3, ...
         (x.^3).*y,  (x.^2).*(y.^2),  x.*(y.^3),  x.^4,  y.^4];
    
  otherwise
    msg = 'ORDER must be 2, 3, or 4.';
    eid = sprintf('Images:%s:invalidOrder',mfilename);
    error(eid,msg);
end

%-------------------------------
% Function  findPiecewiseLinear
%
function [trans,output] = findPiecewiseLinear(uv,xy,options)
%
% For each triangle, finding a plane to warp xy into uv is an affine 
% transformation.
%
% For an affine transformation:
%
%                   [ A D ]
% [u v] = [x y 1] * [ B E ]
%                   [ C F ]
%
% [ u1 v1 ]   [ x1 y1 1 ]   [ A D ]
% [ u2 v2 ] = [ x2 y2 1 ] * [ B E ]
% [ u3 v3 ]   [ x3 y3 1 ]   [ C F ]
%
% Rewriting the above matrix equation:
% U = X * T, where T = [A B C; D E F]'
%
% With the 3 correspondence points of each triangle, we can solve for T,
% T = X\U 
%
% see "Piecewise linear mapping functions for image registration" Ardeshir
% Goshtasby, Pattern Recognition, Vol 19, pp. 459-466, 1986.


% initialize
output.uv = uv;
output.xy = xy;
output.uv_dev = [];
output.xy_dev = [];

x = xy(:,1);
y = xy(:,2);

% Need to pass DELAUNAY options 'QJ' and 'Pp' to avoid qhull errors and 
% warnings for 2-D triangulations.
tri = delaunay(x,y); 

ntri = size(tri,1);

if (ntri<2)
    [msg, eid] = RankError(options.K,'piecewise linear');
    error(eid,msg);
end

% Find all inside-out triangles
bad_triangles =  FindInsideOut(xy,uv,tri);
if ~isempty(bad_triangles)
    
    % find bad_vertices, eliminate bad_vertices
    num_bad_triangles = length(bad_triangles);
    bad_vertices = zeros(num_bad_triangles,1);
    for i = 1:num_bad_triangles
        bad_vertices(i) = FindBadVertex(x,y,tri(bad_triangles(i),:));
    end
    bad_vertices = unique(bad_vertices);
    num_bad_vertices = length(bad_vertices);    

    output.xy_dev = xy(bad_vertices,:); % update to return
    output.uv_dev = uv(bad_vertices,:);
    
    xy(bad_vertices,:) = []; % eliminate bad ones
    uv(bad_vertices,:) = [];
    nvert = size(xy,1);
    
    output.xy = xy; % update to return
    output.uv = uv;

    if (nvert < options.K)
        msg1 = sprintf('Eliminated %d control point pair(s).\n', ...
                       num_bad_vertices);
        msg1 = [msg1 sprintf('Only %d control points remain.\n',nvert)];
        [eid,msg2] = CountError(options.K,'piecewise linear');
        msg = [msg1 msg2];
        error(eid,msg);
    end
    x = xy(:,1);
    y = xy(:,2);
    tri = delaunay(x,y);
    ntri = size(tri,1);
    
    % see if there are any more bad triangles
    more_bad_triangles = FindInsideOut(xy,uv,tri);
    if ~isempty(more_bad_triangles)
        msg = sprintf('Eliminated %d control point pair(s).\n',...
                      num_bad_vertices);
        msg = [msg 'Fold-over triangles remain. See CP2TFORM reference page.'];
        eid = sprintf('Images:%s:foldoverTriangles',mfilename);
        error(eid,msg)
    end
    
    % print report about triangles and how many points were eliminated
    tri_list = sprintf('%d  ', bad_triangles);
    line1 =  ['Fold-over triangle(s):  ' tri_list sprintf('\n')];
    vert_list = sprintf('%d  ', bad_vertices);
    line2 = ['Compensated by eliminating control point pair(s):  ', vert_list];
    warnmsg = [line1 line2];
    wid = 'Images:cp2tform:foldOverTriangles';
    warning(wid, '%s', warnmsg)
    
end

% calculate reverse mapping for each triangle
T = zeros(3,2,ntri);
for itri = 1:ntri
    
    X = [ xy( tri(itri,:), : ) ones(3,1)];
    U =   uv( tri(itri,:), : );
    if (rank(X) >= 3)
        T(:,:,itri) = X\U;
    else
        [msg, eid] = RankError(options.K,'piecewise linear');
        error(eid,msg);
    end

end

% Create TriangleGraph which is a sparse connectivity matrix
nxy = size(xy,1);
S = sparse( repmat((1:ntri)',1,3), tri, 1, ntri, nxy);

% Create OnHull to be 1 for ControlPoints on convex hull, 0 for
% interior points.
hull_indices = convhull(x,y);
OnHull = zeros(size(x));
OnHull(hull_indices) = 1;

tdata.Triangles = tri;
tdata.ControlPoints = xy;
tdata.OnHull = OnHull;
tdata.ConvexHullVertices = hull_indices;
tdata.TriangleGraph = S;
tdata.PiecewiseLinearTData = T;

trans = maketform('custom',2,2,[],@inv_piecewiselinear,tdata);

%-------------------------------
% Function FindInsideOut
%
function index = FindInsideOut(xy,uv,tri)

% look for inside-out triangles using line integrals
x = xy(:,1);
y = xy(:,2);
u = uv(:,1);
v = uv(:,2);

p = size(tri,1);

xx = reshape(x(tri),p,3)';
yy = reshape(y(tri),p,3)';
xysum = sum( (xx([2 3 1],:) - xx).* (yy([2 3 1],:) + yy), 1 );

uu = reshape(u(tri),p,3)';
vv = reshape(v(tri),p,3)';
uvsum = sum( (uu([2 3 1],:) - uu).* (vv([2 3 1],:) + vv), 1 );
    
index = find(xysum.*uvsum<0);

%-------------------------------
% Function FindBadVertex
%
function vertex = FindBadVertex(x,y,vertices)

% Get middle vertex of triangle where "middle" means the largest angle,
% which will have the smallest cosine.

vx = x(vertices)';
vy = y(vertices)';
abc = [ vx - vx([3 1 2]); vy - vy([3 1 2]) ];
a = abc(:,1);
b = abc(:,2);
c = abc(:,3);

% find cosine of angle between 2 vectors
vcos(1) = get_cos(-a, b);
vcos(2) = get_cos(-b, c);
vcos(3) = get_cos( a,-c);

[~, index] = min(vcos);
vertex = vertices(index);

%-------------------------------
% Function get_cos
%
function vcos = get_cos(a,b)

mag_a = sqrt( a(1)*a(1) + a(2)*a(2) );
mag_b = sqrt( b(1)*b(1) + b(2)*b(2) );
vcos = dot(a,b) / (mag_a*mag_b);


%-------------------------------
% Function  findLWM
%
function [trans,output] = findLWM(uv,xy,options)
%
% For a polynomial transformation: 
%
% u = X*A, v = X*B, solve for A and B:
%     A = X\u;
%     B = X\v;   
% 
% The matrix X depends on the order of the polynomial.
% X will be M-by-K, where K = (order+1)*(order+2)/2;
% so A and B will be vectors of length K.
%
%   order = 2   
%     X = [ones(M,1),  x,  y,  x.*y,  x.^2,  y.^2];
%     so X is an M-by-6 matrix 
%
% see "Image registration by local approximation methods" Ardeshir
% Goshtasby, Image and Vision Computing, Vol 6, p. 255-261, 1988.

if (options.order ~= 2)
    msg = 'Internal problem: only polynomials of order=2 are supported.';
    eid = sprintf('Images:%s:internalProblemPolyOrd',mfilename);
    error(eid,msg);
end

output = [];
N = options.N;
M = size(xy,1);

x = xy(:,1);
y = xy(:,2);
u = uv(:,1);
v = uv(:,2);

T = zeros(options.K,2,M);
radii = zeros(M,1);
for icp = 1:M
    
    % find N closest points
    distcp = sqrt( (x-x(icp)).^2 + (y-y(icp)).^2 );
    [dist_sorted,indx] = sort(distcp);
    radii(icp) = dist_sorted(N);
    neighbors = indx(1:N);        
    neighbors = sort(neighbors);
    xcp = x(neighbors);
    ycp = y(neighbors);    
    ucp = u(neighbors);
    vcp = v(neighbors);    

    % set up matrix eqn for polynomial of order=2
    X = [ones(N,1),  xcp,  ycp,  xcp.*ycp,  xcp.^2,  ycp.^2];

    if rank(X)>=options.K
        % u = X*A, v = X*B, solve for A and B:
        A = X\ucp;
        B = X\vcp;
        T(:,:,icp) = [A B];
    else
        [msg, eid] = RankError(options.K,'polynomial');
        error(eid,msg);
    end

end

tdata.LWMTData = T;
tdata.ControlPoints = xy;
tdata.RadiiOfInfluence = radii;

trans = maketform('custom',2,2,[],@inv_lwm,tdata);

%-------------------------------
% Function  ParseInputs
%
function [uv, xy, method, options] = ParseInputs(varargin)

% defaults
options.order = 3;
options.K = [];
N = [];

% iptchecknargin(2,4,nargin,mfilename);

% figure out if syntax is
% CP2TFORM(CPSTRUCT,TRANSFORMTYPE,...) or
% CP2TFORM(INPUT_POINTS,BASE_POINTS,TRANSFORMTYPE,...)

if isa(varargin{1},'struct')
    % TRANS = CP2TFORM(CPSTRUCT,TRANSFORMTYPE)
    % TRANS = CP2TFORM(CPSTRUCT,'polynomial',ORDER)        
    % TRANS = CP2TFORM(CPSTRUCT,'lwm',N)    

    iptchecknargin(2,3,nargin,mfilename);
    
    [uv,xy] = cpstruct2pairs(varargin{1});
    method = getMethod(varargin{2});

    nargs_to_go = nargin - 2;
    if nargs_to_go > 0
        args = varargin(3:end);
    end

else
    % TRANS = CP2TFORM(INPUT_POINTS,BASE_POINTS,TRANSFORMTYPE)
    % TRANS = CP2TFORM(INPUT_POINTS,BASE_POINTS,'polynomial',ORDER)        
    % TRANS = CP2TFORM(INPUT_POINTS,BASE_POINTS,'lwm',N)    

%     iptchecknargin(3,4,nargin,mfilename);
    
    uv = varargin{1};
    xy = varargin{2};
    method = getMethod(varargin{3});

    nargs_to_go = nargin - 3;
    if nargs_to_go > 0
        args = varargin(4:end);
    end
    
end
    
if size(uv,2) ~= 2 || size(xy,2) ~= 2
    eid = sprintf('Images:%s:invalidControlPointMatrix', mfilename);
    error(eid,...
          'In function %s, control point matrices must be M-by-2.', mfilename);
end

if size(uv,1) ~= size(xy,1)
    eid = sprintf('Images:%s:needSameNumControlPoints', mfilename);
    error(eid,...
          'In function %s, INPUT and BASE images need same number of control points.', mfilename);
end

switch nargs_to_go
  case 0
    % TRANS = CP2TFORM(INPUT_POINTS,BASE_POINTS,TRANSFORMTYPE)
    % TRANS = CP2TFORM(CPSTRUCT,TRANSFORMTYPE)
    
  case 1
    if strcmp(method,'polynomial')
        % TRANS = CP2TFORM(INPUT_POINTS,BASE_POINTS,'polynomial',ORDER)        
        % TRANS = CP2TFORM(CPSTRUCT,'polynomial',ORDER)
        options.order = args{1};    

    elseif strcmp(method,'lwm')
        % TRANS = CP2TFORM(INPUT_POINTS,BASE_POINTS,'lwm',N)        
        % TRANS = CP2TFORM(CPSTRUCT,'lwm',N)                
        N = args{1};    

    else 
        eid = sprintf('Images:%s:tooManyInputs', mfilename);
        error(eid,...
              'Too many input arguments were passed to %s.', mfilename);
    end
    
  otherwise
    eid = sprintf('Images:%s:tooManyInputs', mfilename);
    error(eid,...
          'Too many input arguments were passed to %s.', mfilename);
end

switch method
  case 'nonreflective similarity'
    options.K = 2;
  case 'similarity'
    options.K = 3; 
  case 'affine'
    options.K = 3; 
  case 'projective'
    options.K = 4;
  case 'polynomial'
    order = options.order;
    % validate order
    if ~isnumeric(order) || numel(order)~=1 || sum(order==2:4)~=1
        eid = sprintf('Images:%s:invalidPolynomialOrder', mfilename);
        error(eid,...
              'In function %s, polynomial order must be 2, 3, 4.', mfilename);
    end
    options.K = (order+1)*(order+2)/2;
  case 'piecewise linear'
    options.K = 4; % for 'piecewise linear' need at least 4 points, 2 triangles
  case 'lwm'
    order = 2;
    options.K = (order+1)*(order+2)/2;
    options.order = order;

    if isempty(N)
        % conservative default N protects user from ill-conditioned polynomials
        N = 2*options.K; 
                         
    else
        % validate N
        if ~isnumeric(N) || numel(N)~=1 || rem(N,1)~=0 || N<options.K
            eid = sprintf('Images:%s:invalidInputN', mfilename);
            error(eid,...
                  'In function %s, N must be an integer greater than or equal to %d.',...
                          mfilename, options.K);
        end
    end
    options.N = N;
  otherwise
    eid = sprintf('Images:%s:internalProblem', mfilename);
    error(eid,...
          'Function %s has an internal problem: unrecognized method.',...
          mfilename);
    
end    

%-------------------------------
% Function  getMethod
%
function method = getMethod(method_string)

    method_string = lower(method_string);
    
    % Figure out which method to use
    methods = {'affine',...
               'nonreflective similarity',...
               'linear conformal',... % synonym for 'nonreflective similarity'
               'similarity',...
               'projective', 'polynomial', ...
               'piecewise linear', 'lwm'};
    if ischar(method_string)
        indx = strmatch(method_string, methods);
        switch length(indx)
          case 0
            eid = sprintf('Images:%s:unrecognizedTransformType',mfilename);
            error(eid,'Unrecognized TRANSFORMTYPE ''%s'' in function %s.',method_string,mfilename);
          case 1
            method = methods{indx};
          otherwise
            eid = sprintf('Images:%s:ambiguousTransformType',mfilename);
            error(eid,'Ambiguous TRANSFORMTYPE ''%s'' in function %s.',method_string,mfilename);
        end
    else
        eid = sprintf('Images:%s:transformTypeIsNotString',mfilename);
        error(eid,'TRANSFORMTYPE must be a string in function %s.',mfilename);
    end	

    if strcmp(method,'linear conformal')
        method = 'nonreflective similarity';
    end
        
%-------------------------------
% Function  CountError
%
function [msg, eid] = CountError(K,transform_string)

msg = sprintf('At least %d points needed to infer %s transform.', ...
              K,transform_string);
eid = sprintf('Images:%s:atLeast%dPointsReq',mfilename,K);

%-------------------------------
% Function  RankError
%
function [msg, eid] = RankError(K,transform_string)

msg = 'At least %d non-collinear points needed to infer %s transform.';
msg = sprintf(msg,K,transform_string);
eid = sprintf('Images:%s:atLeast%dNonColinearPointsReq',mfilename,K);
