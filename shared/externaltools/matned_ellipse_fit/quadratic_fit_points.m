function quad = quadratic_fit_points(pts, iso, make_plot)

% quad = quadratic_fit_points(pts, iso)
%
% Input:
%  'pts' should be a Nx3 matrix where each line contains an entry (x, y, w),
%   where 'w' is the weight associated with the position (x,y). If all the
%   weights are equal, it is sufficient to provide a Nx2 matrix of (x, y)
%
% Output:
%  coefficients q = [ a, b, c, d, e, f ] of the functional
%  quad(X,Y) = q(1)*X^2 + q(2)*X*Y + q(3)*Y^2 + q(4)*X + q(5)*Y + q(6)
%  that is minimal the summed squared error at the points.
%
%
% if ( iso == 1 ) the fit is constrained to be circular, with a==b.
%
% The characteristics of the ellipse corresponding to the coefficients in 'quad'
% can be calculated by 'quadratic_center', and the ellipse can be plotted
% on top of an existing plot with 'quadratic_plot'
%
% F. Nedelec, April 2014

if nargin < 3
    make_plot = 0;
end

if nargin < 2
    iso = 0;
end

%%

x = double(pts(:,1));
y = double(pts(:,2));

xx = x .* x;
xy = x .* y;
yy = y .* y;

%% various sums:

if size(pts, 2) == 2
    
    sw     = size(x, 1);
    sx     = sum( x );
    sy     = sum( y );
    
    sxx    = sum( xx );
    sxy    = sum( xy );
    syy    = sum( yy );
    
    sxxx   = sum( xx .* x );
    sxxy   = sum( xx .* y );
    sxyy   = sum( xy .* y );
    syyy   = sum( yy .* y );
    
    sxxxx  = sum( xx .* xx );
    sxxxy  = sum( xx .* xy );
    sxxyy  = sum( xx .* yy );
    sxyyy  = sum( xy .* yy );
    syyyy  = sum( yy .* yy );
    
else
    
    w = double(pts(:,3));
    
    sw     = sum( w );
    
    if sw <= 0
        error('the sum of the weights is negative');
    end
    
    sx     = sum( x .* w );
    sy     = sum( y .* w );
    
    sxx    = sum( xx .* w );
    sxy    = sum( xy .* w );
    syy    = sum( yy .* w );
    
    sxxx   = sum( xx .* x .* w );
    sxxy   = sum( xx .* y .* w );
    sxyy   = sum( xy .* y .* w );
    syyy   = sum( yy .* y .* w );
    
    sxxxx  = sum( xx .* xx .* w );
    sxxxy  = sum( xx .* xy .* w );
    sxxyy  = sum( xx .* yy .* w );
    sxyyy  = sum( xy .* yy .* w );
    syyyy  = sum( yy .* yy .* w );

end

%% matrix 

if iso > 0
   
    S = [ sxx, sxy, sx; ...
          sxy, syy, sy; ...
          sx,  sy,  sw ];

    c = S \ [ -sxxx-sxyy; -syyy-sxxy; -sxx-syy ];

    quad = [ 1, 0, 1, c(1), c(2), c(3) ];

else
    
    if ( 0 )
        % this assumes that you can normalize the equation by f=1,
        % it fails if the solution is a circle that passes through the origin
        S = [ sxxxx, sxxxy, sxxyy, sxxx, sxxy; ...
              sxxxy, sxxyy, sxyyy, sxxy, sxyy; ...
              sxxyy, sxyyy, syyyy, sxyy, syyy; ...
              sxxx,  sxxy,  sxyy,  sxx,  sxy ; ...
              sxxy,  sxyy,  syyy,  sxy,  syy ];

        c = S \ [ sxx; sxy; syy; sx; sy ];
        quad = [ c(1), c(2), c(3), c(4), c(5), -1 ];
    end
    
    % the matrix is singular, and the nullspace vector gives the solution
    M = [ sxxxx, sxxxy, sxxyy, sxxx, sxxy, sxx; ...
          sxxxy, sxxyy, sxyyy, sxxy, sxyy, sxy; ...
          sxxyy, sxyyy, syyyy, sxyy, syyy, syy; ...
          sxxx,  sxxy,  sxyy,  sxx,  sxy,  sx;  ...
          sxxy,  sxyy,  syyy,  sxy,  syy,  sy;  ...
          sxx,   sxy,   syy,   sx,   sy,   sw];

    [ eigvec, ~ ] = eig(M);
    
    quad = eigvec(:,1)';
    quad = quad .* sqrt( 2 /( quad(1)^2 + quad(3)^2 ));
    
end
    
%%

%fprintf(' Quadratic fit: %.3f X^2  %+.3f XY  %+.3f Y^2  %+.3f X  %+.3f Y  %+.3f\n', quad);
    
if make_plot > 0
    
    figure;
    hold on;
    if exist('w', 'var') == 1
        m = max(w);
        col = w/m;
    else
        col = 0.5 * ones(size(pts,1), 3);
    end
        
    for n = 1:length(pts)
        cx = pts(n,1);
        cy = pts(n,2);
        r = 0.5;
        patch([ cx-r, cx+r, cx+r, cx-r ], [ cy+r, cy+r, cy-r, cy-r ], col(n), 'EdgeColor', 'none');
    end
    colormap('gray');
    axis equal;
    quadratic_plot(quad);
    
end

end

