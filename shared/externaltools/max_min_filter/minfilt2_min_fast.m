function Y = minfilt2_min_fast(X,S)
%  MINFILT2    Two-dimensional min filter
%
%     Y = MINFILT2(X,[M N]) performs two-dimensional minimum
%     filtering on the image X using an M-by-N window. The result
%     Y contains the minimun value in the M-by-N neighborhood around
%     each pixel in the original image. 
%     This function uses the van Herk algorithm for min filters.
%
%     Y = MINFILT2(X,M) is the same as Y = MINFILT2(X,[M M])
%
%     Y = MINFILT2(X) uses a 3-by-3 neighborhood.
%
%     Y = MINFILT2(..., 'shape') returns a subsection of the 2D
%     filtering specified by 'shape' :
%        'full'  - Returns the full filtering result,
%        'same'  - (default) Returns the central filter area that is the
%                   same size as X,
%        'valid' - Returns only the area where no filter elements are outside
%                  the image.
%
%     See also : MAXFILT2, VANHERK
%

% Initialization
% [S, shape] = parse_inputs(varargin{:});

if nargin<2
    S=[3 3];
elseif length(S) == 1
   S(2) = S(1);
end
% filtering
Y = vanherk_min_fast(X,S(1),'row');
Y = vanherk_min_fast(Y,S(2),'col');
if size(Y,1)<size(X,1)
    Y(end+1,:)=0;
end
if size(Y,2)<size(X,2)
    Y(:,end+1)=0;
end

