function quad = fit_ellipse(im, threshold, make_plot)

% function quad = fit_ellipse(image, threshold, make_plot)
% function quad = fit_ellipse(movie, threshold, make_plot)
%
% Fit a quadratic form to the part of the image that is above the threshold
% or to each image in the movie. 
%
% F. Nedelec, June 2014

if nargin < 3
    make_plot = 0;
end

%%


if isfield(im, 'data')

    quad = zeros(length(im), 6);
    for i = 1:length(im)
        quad(i, :) = quadratic_fit_points(make_points(im(i).data, threshold), 0, make_plot-1);
    end
    
else

    % here '0' specifies non-isotropic
    quad = quadratic_fit_points(make_points(im, threshold), 0, make_plot-1);
    
    if make_plot > 0
        show_image(im);
        quadratic_plot(quad, 1);
    end
    
end


end

    function pts = make_points(im, threshold)
        [x, y] = find( im > threshold );
        ind = logical( im > threshold );
        w = im(ind) - threshold;
        pts = cat(2, x, y, w);
    end
