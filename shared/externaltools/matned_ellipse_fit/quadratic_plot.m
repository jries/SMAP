function h = quadratic_plot(quad, flip)

% function h = quadratic_plot(q, flip)
%
% Draw the ellipse specified by the argument:
% q(1)*X^2 + q(2)*X*Y + q(3)*Y^2 + q(4)*X + q(5)*Y + q(6) = 0,
%  
%
% if flip=1, the X and Y coordinates are flipped, which is necessary to
% overlay on top of images.
%
% F. Nedelec, April 2014


if nargin < 2
    flip = 0;
end

[cen, axs, angle] = quadratic_center(quad);

cx = cen(1);
cy = cen(2);

ca = cos(angle);
sa = sin(angle);

theta = 0 : 0.1 : 2*pi+0.1;

px = cx + ca*axs(1)*cos(theta) - sa*axs(2)*sin(theta);
py = cy + sa*axs(1)*cos(theta) + ca*axs(2)*sin(theta);

if ( flip )
    h = plot(py, px);
    hold on;
    plot(cy, cy, 'o');
    plot([cy+sa*axs(1), cy-sa*axs(1)], [cx+ca*axs(1), cx-ca*axs(1)]);
else
    h = plot(px, py);
    hold on;
    plot(cx, cy, 'o');
    plot([cx+ca*axs(1), cx-ca*axs(1)], [cy+sa*axs(1), cy-sa*axs(1)]);
end

end

