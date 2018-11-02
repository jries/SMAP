function [ Voronoi ] = VorArea( A1, ko, ke )
% Performs Voronoi tessellation on super-res data and calculates the 
% inverted areas for Voronoi cells
% Optional parameters: fov, ko, ke

if isstruct(A1)
A1 = A1.data;
end
if ~exist('ko', 'var')
    ko = 1;
end
if ~exist('ke', 'var')
    ke = size(A1, 1);
end

A1 = round(A1(ko:ke, :));
Anew = A1(:, 4:5);
[A, ~, ic] = unique(Anew,'rows');
[V,C] = voronoin(A);
Area = zeros(size(C));
for i = 1 : size(C, 1)
    Area(i) = polyarea(V(C{i}, 1), V(C{i}, 2));
end
Area = 1./Area;
[n, ~] = histc(ic, unique(ic)); % n = number of repeats for each point with index ic.
Area = Area .* n;
Area(isnan(Area) | isinf(Area)) = 0;
Voronoi = { Area, V, C, A1, ic };
