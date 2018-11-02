function handle = draw_ellipse(cen, ell, hAxes)

% handle = draw_ellipse(cen, ell)
% handle = draw_ellipse(cen, ell, axes_handle)
%
% Draw an ellipse in the current graph, centered on cen,
% and defined by ell:
%   ell(1) = angle of X with major axis,
%   ell(2) = half-length of major axis
%   ell(3) = half-length of minor axis
%
% returns a handle to the graphic object created
%
% F. Nedelec, may 2000; Email: nedelec@embl.de


if length(cen) ~= 2  ||  ~isreal(cen)
    disp(cen);
    error('First argument should be a real vector of length 2');
end

if length(ell) ~= 3  ||  ~isreal(ell)
    disp(ell);
    error('Second argument should be a real vector of length 3');
end

if all( ell == 0 )
    error('draw_ellipse:null', 'Null ellipse arguments');
end


tilt = ell(1);
rotm = [ cos(tilt), -sin(tilt); sin(tilt), cos(tilt) ];

theta = 0:0.1:2*pi;
pts = cen' * ones(size(an)) + rotm * [ ell(2)*cos(theta); ell(3)*sin(theta) ];

if nargin == 3
    handle = plot(hAxes, pts(2,:), pts(1,:), 'y.', 'MarkerSize', 1);
else
    handle = plot(pts(2,:), pts(1,:), 'y.', 'MarkerSize', 1);
end

end