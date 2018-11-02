%ISECT   Find intersections of functions a and b defined at gridpoints x
%
% SYNOPSIS:
%   intersections = isect(x,a,b)

% (C) Copyright 2012               Quantitative Imaging Group
%     All rights reserved          Faculty of Applied Physics
%                                  Delft University of Technology
%                                  Lorentzweg 1
%                                  2628 CJ Delft
%                                  The Netherlands
% Robert Nieuwenhuizen, Oct 2012

function out = isect(x,a,b)

if (size(x)~= size(a)) OR (size(x)~= size(b))
    error('x, a, and b must have the same size.');
end

% Check for uniqueness of x-values
if size(x)~= size(unique(x))
    error('x must have unique values.');
end

if isrow(x)
    x = x';
    a = a';
    b = b';
end

% Remove NaN points from arrays
sortmat = [x a b];
sortmat = sortmat(~isnan(a),:);
sortmat = sortmat(~isnan(b),:);

% Sort a,b and x based on x
sortmat = sortrows(sortmat);
x = sortmat(:,1);
a = sortmat(:,2);
b = sortmat(:,3);
clear sortmat

% Calculate the difference between a and b
a = a-b;

% Indices for the values of x before a and b intersect
i_inds = ((a(2:end)>0)&(a(1:(end-1))<0))|((a(2:end)<0)&(a(1:(end-1))>0));

% x-values for points not in x where a and b intersect 
i_points = i_inds.*(x(1:(end-1)) - a(1:(end-1))./(a(2:end)-a(1:(end-1))).*(x(2:end)-x(1:(end-1))));

% Include gridpoints where a and b intersect
if sum(a==0)~=0
    out = cat(1,x(a==0),i_points(i_inds));
else
    out = i_points(i_inds);
end

% Remove undesirable values
out = out(~isnan(out));
out = out(~isinf(out));

% Sort output
out = sort(out);

end
