function [cens, axes] = ellipse_radius(mov, level)

% Extract center and axes of radius for each image in a movie
%
% function ellipse_radius(mov, level)
%
% F. Nedelec, June 2014

nbi = length(mov);

%% 

if nargin < 2
    
    s = double(mov(1).data);
    for i = 2:nbi
        s = s + double(mov(i).data);
    end
    s = s / nbi;
    [level, sig] = image_background(s);
    level = level + 2 * sig;
    level = manual_threshold(s, level);
    fprintf('Selected threshold = %f\n', level);

end

%% 

quads = fit_ellipse(mov, level);

cens = [];
axes = [];

for n = 1:nbi
    [cen, axs, angle] = quadratic_center(quads(n,:));
    cens = cat(1, cens, cen);
    axes = cat(1, axes, axs);
end

end

