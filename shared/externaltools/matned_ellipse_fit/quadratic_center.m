function [ cen, axs, angle ] = quadratic_center(q, verbose)

% function [ cen, axes, angle ] = quadratic_center( q )
%
% returns the center, axes length and tilt of the ellipse defined by the equation:
%  q(1)*X^2 + q(2)*X*Y + q(3)*Y^2 + q(4)*X + q(5)*Y + q(6) = 0
%
% F. Nedelec, April 2014

if length(q) ~= 6
    error('input must be a vector of size 6');
end

if nargin < 2
    verbose = 0;
end

if verbose
    fprintf(' Quadratic    : %.3f X^2  %+.3f XY  %+.3f Y^2  %+.3f X  %+.3f Y  %+.3f\n', q);
end

%% Center

a = q(1);
b = q(2) / 2;
c = q(3);
d = q(4) / 2;
f = q(5) / 2;
g = q(6);

det = b*b - a*c;

if det > 0
    error('this is not an ellipse');
end

cen = [ c*d-b*f,  a*f-b*d ] ./ det;

%% rotate formula to align with axes

if nargout > 1
    
    [ qr, angle ] = reduce_quadratic(q);
    % center of rotated ellipse:
    cx = -0.5*qr(4)/qr(1);
    cy = -0.5*qr(5)/qr(3);

    % calculate first axis with Y = cy:
    yy = ( qr(3)*cy + qr(5) ) * cy + qr(6);
    ax1 =  0.5 * sqrt( qr(4).^2 - 4*qr(1)*yy ) / abs(qr(1));
    
    % calculate first axis with X = cx:
    xx = ( qr(1)*cx + qr(4) ) * cx + qr(6);
    ax2 =  0.5 * sqrt( qr(5).^2 - 4*qr(3)*xx ) / abs(qr(3));
    
    if ( ax2 > ax1 )
        axs = [ ax2, ax1 ];
        angle = angle + pi/2;
    else
        axs = [ ax1, ax2 ];
    end
    

    if verbose
        fprintf(' Reduced form : %.3f X^2  %+.3f XY  %+.3f Y^2  %+.3f X  %+.3f Y  %+.3f\n', qr);
        fprintf(' Angle        : %.3f\n', angle);
        fprintf(' Center       : %.3f %.3f\n', cx, cy);
        fprintf(' Axes         : %.3f %.3f\n', ax1, ax2);
    end
end

return
%% Alternative method 

if nargout > 1

    nu = a*f*f + c*d*d + g*b*b - 2*b*d*f - a*c*g;
    dn = sqrt( (a-c)^2 + 4*b*b );
    ax1 = sqrt( 2 * nu / ( det * (  dn - a - c ) ));
    ax2 = sqrt( 2 * nu / ( det * (- dn - a - c ) ));
    axl = [ ax1, ax2 ];

end

if nargout > 2
    
    if b == 0
        if a <= c
            angle = 0;
        else
            angle = pi/2;
        end
    else
        if a < c
            angle = 0.5 * atan2(2*b, a-c);
        else
            angle = 0.5 * ( pi + atan2(2*b, a-c) );
        end
    end 
    
end

end

%%

function [ qr, angle ] = reduce_quadratic(q)

    angle = 0.5 * atan(q(2)/(q(1)-q(3)));
    
    if isnan(angle)
        angle = 0;
    end
    
    qr = rotate_quadratic(q, angle);

end


function qr = rotate_quadratic(q, angle)

    % Apply a rotation of angle a:
    %   X = cos(a)*x - sin(a)*y
    %   Y = sin(a)*x + cos(a)*y
    ca = cos(angle);
    sa = sin(angle);

    qr(1) = q(1)*ca*ca + q(2)*ca*sa + q(3)*sa*sa;
    qr(2) = q(2)*(ca*ca-sa*sa) + 2*ca*sa*(q(3)-q(1));
    qr(3) = q(1)*sa*sa - q(2)*sa*ca + q(3)*ca*ca;
    qr(4) = q(4)*ca + q(5)*sa;
    qr(5) = q(5)*ca - q(4)*sa;
    qr(6) = q(6);
    
end

