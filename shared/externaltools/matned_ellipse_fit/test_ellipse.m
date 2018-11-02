% a script to test the ellipse with one of Aastha's picture
% F. Nedelec, 15 October 2015

im = imread('platelet.tif');

% extract one channel:
i = im(:,:,1);

% threshold value is 700 (this parameter is critical)
fit_ellipse(i, 700, 1);

