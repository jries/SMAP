function [ I, Avor ] = drawVor( A1, p, fov, ko, ke )
% Builds the super-res image with help of Voronoi triangulation. The
% density at each data point is estimated as an inverse of the area of the
% corresponding Voronoi cell. This value is added to the 10-th column of
% the original eventlist yelding 'Avor'. 
% The values of the density are then interpolated to an regular grid with a
% step size of p (nm) using the method 'natural'.
% Optional parameters: fov, ko, ke
if ~exist('ko', 'var')
    ko = 1;
end
if ~exist('ke', 'var')
    ke = size(A1, 1);
end
Vor = VorArea( A1, ko, ke );
Area = Vor{1};
A1 = Vor{4};
ic = Vor{5};
Areacorr = zeros(size(A1,1),1);
Areacorr(:) = Area(ic,:);
Avor = [A1, Areacorr];
%Avor(:,10) = Avor(:,10) .* Avor(:,7) ./ mean(Avor(:,7)); %multiply by number of photons
% drawing
if exist ('p', 'var')
s = round(fov / p); % size in pixels
M = zeros(s * s, 3); %template for query matrix
a = 1;
for i = 1:s
    for j = 1:s
        M(a,1) = j;
        M(a,2) = i;
        a = a+1;
    end
end
M = M * p; %convert to nanometers

% interpolate Anew to M - v4 is too slow, natural OK
M(:,3) = griddata(Avor(:,4), Avor(:,5), Avor(:,10), M(:,1), M(:,2), 'natural');

%delete negative values
M(M(:,3) < 0 , 3) = 0;
M(isnan(M(:,3)), 3) = 0;
M(isinf(M(:,3)), 3) = 0;

% build an image from M(x,y,z)
I = zeros(s);
for k = 1 : size(M,1)
%     I(M(k,2) / p, M(k,1) / p) = M(k,3);
    I(round(M(k,2) / p), round(M(k,1) / p)) = M(k,3);
end
else
    I = [];
end
end
%figure; imshow(10*uint8(256*I/max(max(I)))); colormap(hot);

