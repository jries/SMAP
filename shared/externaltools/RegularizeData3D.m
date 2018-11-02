function [zgrid,xgrid,ygrid] = RegularizeData3D(x,y,z,xnodes,ynodes,varargin)
% RegularizeData3D: Produces a smooth 3D surface from scattered input data.
%
%          RegularizeData3D is a modified version of GridFit from the Matlab File Exchange.
%          RegularizeData3D does essentially the same thing, but is an attempt to overcome several
%          shortcomings inherent in the design of the legacy code in GridFit.
%
%          * GridFit lacks cubic interpolation capabilities.
%            Interpolation is necessary to map the scattered input data to locations on the output
%            surface.  The output surface is most likely nonlinear, so linear interpolation used in
%            GridFit is a lousy approximation.  Cubic interpolation accounts for surface curvature,
%            which is especially beneficial when the output grid is coarse in x and y.
%
%          * GridFit's "smoothness" parameter was poorly defined and its methodology may have led to bad output data.
%            In RegularizeData3D the smoothness parameter is actually the ratio of smoothness (flatness) to
%            fidelity (goodness of fit) and is not affected by the resolution of the output grid.
%            Smoothness = 100 gives 100 times as much weight to smoothness (and produces a nearly flat output
%                         surface)
%            Smoothness = 1 gives equal weight to smoothness and fidelity (and results in noticeable smoothing)
%						 Smoothness = 0.01 gives 100 times as much weight to fitting the surface to the scattered input data (and
%						              results in very little smoothing)
%						 Smoothness = 0.001 is good for data with low noise.  The input points nearly coincide with the output
%                         surface.
%
%          * GridFit didn't do a good job explaining what math it was doing; it just gave usage examples.
%            
%            For a detailed explanation of "the math behind the magic" on a 3D dataset, see:
%            http://mathformeremortals.wordpress.com/2013/07/22/introduction-to-regularizing-3d-data-part-1-of-2/
%
%            and to apply the same principles to a 2D dataset see:
%            http://mathformeremortals.wordpress.com/2013/01/29/introduction-to-regularizing-with-2d-data-part-1-of-3/
%
%            Both of these links include Excel spreadsheets that break down the calculations step by step
%            so that you can see how RegularizeData3D works.  There are also very limited (and very slow!) Excel
%            spreadsheet functions that do the same thing in 2D or 3D.
%
%
%          Aside from the above changes, most of the GridFit code is left intact.
%          The original GridFit page is:
%          http://www.mathworks.com/matlabcentral/fileexchange/8998-surface-fitting-using-gridfit
% usage #1: zgrid = RegularizeData3D(x, y, z, xnodes, ynodes);
% usage #2: [zgrid, xgrid, ygrid] = RegularizeData3D(x, y, z, xnodes, ynodes);
% usage #3: zgrid = RegularizeData3D(x, y, z, xnodes, ynodes, prop, val, prop, val,...);
%
% Arguments: (input)
%  x,y,z - vectors of equal lengths, containing arbitrary scattered data
%          The only constraint on x and y is they cannot ALL fall on a
%          single line in the x-y plane. Replicate points will be treated
%          in a least squares sense.
%
%          ANY points containing a NaN are ignored in the estimation
%
%  xnodes - vector defining the nodes in the grid in the independent
%          variable (x). xnodes need not be equally spaced. xnodes
%          must completely span the data. If they do not, then the
%          'extend' property is applied, adjusting the first and last
%          nodes to be extended as necessary. See below for a complete
%          description of the 'extend' property.
%
%          If xnodes is a scalar integer, then it specifies the number
%          of equally spaced nodes between the min and max of the data.
%
%  ynodes - vector defining the nodes in the grid in the independent
%          variable (y). ynodes need not be equally spaced.
%
%          If ynodes is a scalar integer, then it specifies the number
%          of equally spaced nodes between the min and max of the data.
%
%          Also see the extend property.
%
%  Additional arguments follow in the form of property/value pairs.
%  Valid properties are:
%    'smoothness', 'interp', 'solver', 'maxiter'
%    'extend', 'tilesize', 'overlap'
%
%  Any UNAMBIGUOUS shortening (even down to a single letter) is
%  valid for property names. All properties have default values,
%  chosen (I hope) to give a reasonable result out of the box.
%
%   'smoothness' - scalar or vector of length 2 - the ratio of
%          smoothness to fidelity of the output surface. This must be a
%          positive real number.
%
%          A smoothness of 1 gives equal weight to fidelity (goodness of fit)
%          and smoothness of the output surface.  This results in noticeable
%          smoothing.  If your input data x,y,z have little or no noise, use
%          0.01 to give smoothness 1% as much weight as goodness of fit.
%          0.1 applies a little bit of smoothing to the output surface.
%
%          If this parameter is a vector of length 2, then it defines
%          the relative smoothing to be associated with the x and y
%          variables. This allows the user to apply a different amount
%          of smoothing in the x dimension compared to the y dimension.
%
%          DEFAULT: 0.01
%
%
%   'interp' - character, denotes the interpolation scheme used
%          to interpolate the data.
%
%          DEFAULT: 'triangle'
%
%          'bicubic' - use bicubic interpolation within the grid
%                     This is the most accurate because it accounts
%                     for the fact that the output surface is not flat.
%                     In some cases it may be slower than the other methods.
%
%          'bilinear' - use bilinear interpolation within the grid
%
%          'triangle' - split each cell in the grid into a triangle,
%                      then apply linear interpolation inside each triangle
%
%          'nearest' - nearest neighbor interpolation. This will
%                     rarely be a good choice, but I included it
%                     as an option for completeness.
%
%
%   'solver' - character flag - denotes the solver used for the
%          resulting linear system. Different solvers will have
%          different solution times depending upon the specific
%          problem to be solved. Up to a certain size grid, the
%          direct \ solver will often be speedy, until memory
%          swaps causes problems.
%
%          What solver should you use? Problems with a significant
%          amount of extrapolation should avoid lsqr. \ may be
%          best numerically for small smoothnesss parameters and
%          high extents of extrapolation.
%
%          Large numbers of points will slow down the direct
%          \, but when applied to the normal equations, \ can be
%          quite fast. Since the equations generated by these
%          methods will tend to be well conditioned, the normal
%          equations are not a bad choice of method to use. Beware
%          when a small smoothing parameter is used, since this will
%          make the equations less well conditioned.
%
%          DEFAULT: 'normal'
%
%          '\' - uses matlab's backslash operator to solve the sparse
%                     system. 'backslash' is an alternate name.
%
%          'symmlq' - uses matlab's iterative symmlq solver
%
%          'lsqr' - uses matlab's iterative lsqr solver
%
%          'normal' - uses \ to solve the normal equations.
%
%
%   'maxiter' - only applies to iterative solvers - defines the
%          maximum number of iterations for an iterative solver
%
%          DEFAULT: min(10000,length(xnodes)*length(ynodes))
%
%
%   'extend' - character flag - controls whether the first and last
%          nodes in each dimension are allowed to be adjusted to
%          bound the data, and whether the user will be warned if
%          this was deemed necessary to happen.
%
%          DEFAULT: 'warning'
%
%          'warning' - Adjust the first and/or last node in
%                     x or y if the nodes do not FULLY contain
%                     the data. Issue a warning message to this
%                     effect, telling the amount of adjustment
%                     applied.
%
%          'never'  - Issue an error message when the nodes do
%                     not absolutely contain the data.
%
%          'always' - automatically adjust the first and last
%                     nodes in each dimension if necessary.
%                     No warning is given when this option is set.
%
%
%   'tilesize' - grids which are simply too large to solve for
%          in one single estimation step can be built as a set
%          of tiles. For example, a 1000x1000 grid will require
%          the estimation of 1e6 unknowns. This is likely to
%          require more memory (and time) than you have available.
%          But if your data is dense enough, then you can model
%          it locally using smaller tiles of the grid.
%
%          My recommendation for a reasonable tilesize is
%          roughly 100 to 200. Tiles of this size take only
%          a few seconds to solve normally, so the entire grid
%          can be modeled in a finite amount of time. The minimum
%          tilesize can never be less than 3, although even this
%          size tile is so small as to be ridiculous.
%
%          If your data is so sparse than some tiles contain
%          insufficient data to model, then those tiles will
%          be left as NaNs.
%
%          DEFAULT: inf
%
%
%   'overlap' - Tiles in a grid have some overlap, so they
%          can minimize any problems along the edge of a tile.
%          In this overlapped region, the grid is built using a
%          bi-linear combination of the overlapping tiles.
%
%          The overlap is specified as a fraction of the tile
%          size, so an overlap of 0.20 means there will be a 20%
%          overlap of successive tiles. I do allow a zero overlap,
%          but it must be no more than 1/2.
%
%          0 <= overlap <= 0.5
%
%          Overlap is ignored if the tilesize is greater than the
%          number of nodes in both directions.
%
%          DEFAULT: 0.20
%
%
% Arguments: (output)
%  zgrid   - (nx,ny) array containing the fitted surface
%
%  xgrid, ygrid - as returned by meshgrid(xnodes,ynodes)
%
%
% Speed considerations:
%  Remember that gridfit must solve a LARGE system of linear
%  equations. There will be as many unknowns as the total
%  number of nodes in the final lattice. While these equations
%  may be sparse, solving a system of 10000 equations may take
%  a second or so. Very large problems may benefit from the
%  iterative solvers or from tiling.
%
%
% Example usage:
%
%  x = rand(100,1);
%  y = rand(100,1);
%  z = exp(x+2*y);
%  xnodes = 0:.1:1;
%  ynodes = 0:.1:1;
%
%  g = RegularizeData3D(x,y,z,xnodes,ynodes);
%
% Note: this is equivalent to the following call:
%
%  g = RegularizeData3D(x,y,z,xnodes,ynodes, ...
%              'smooth',1, ...
%              'interp','triangle', ...
%              'solver','normal', ...
%              'gradient', ...
%              'extend','warning', ...
%              'tilesize',inf);
%
%
% Rereleased with improvements as RegularizeData3D
%   2014
%   - Added bicubic interpolation
%   - Fixed a bug that caused smoothness to depend on grid fidelity
%   - Removed the "regularizer" setting and documented the calculation process
% Original Version:
%   Author: John D'Errico
%   e-mail address: woodchips@rochester.rr.com
%   Release: 2.0
%   Release date: 5/23/06

% set defaults
% The default smoothness is 0.01.  i.e. assume the input data x,y,z
% have little or no noise.  This is different from the legacy code,
% which used a default of 1.
params.smoothness = 0.01;
params.interp = 'triangle';
params.solver = 'backslash';
params.maxiter = [];
params.extend = 'warning';
params.tilesize = inf;
params.overlap = 0.20;
params.mask = []; 

% was the params struct supplied?
if ~isempty(varargin)
  if isstruct(varargin{1})
    % params is only supplied if its a call from tiled_gridfit
    params = varargin{1};
    if length(varargin)>1
      % check for any overrides
      params = parse_pv_pairs(params,varargin{2:end});
    end
  else
    % check for any overrides of the defaults
    params = parse_pv_pairs(params,varargin);

  end
end

% check the parameters for acceptability
params = check_params(params);

% ensure all of x,y,z,xnodes,ynodes are column vectors,
% also drop any NaN data
x=x(:);
y=y(:);
z=z(:);
k = isnan(x) | isnan(y) | isnan(z);
if any(k)
  x(k)=[];
  y(k)=[];
  z(k)=[];
end
xmin = min(x);
xmax = max(x);
ymin = min(y);
ymax = max(y);

% did they supply a scalar for the nodes?
if length(xnodes)==1
  xnodes = linspace(xmin,xmax,xnodes)';
  xnodes(end) = xmax; % make sure it hits the max
end
if length(ynodes)==1
  ynodes = linspace(ymin,ymax,ynodes)';
  ynodes(end) = ymax; % make sure it hits the max
end

xnodes=xnodes(:);
ynodes=ynodes(:);
dx = diff(xnodes);
dy = diff(ynodes);
nx = length(xnodes);
ny = length(ynodes);
ngrid = nx*ny;

% check to see if any tiling is necessary
if (params.tilesize < max(nx,ny))
  % split it into smaller tiles. compute zgrid and ygrid
  % at the very end if requested
  zgrid = tiled_gridfit(x,y,z,xnodes,ynodes,params);
else
  % its a single tile.
  
  % mask must be either an empty array, or a boolean
  % aray of the same size as the final grid.
  nmask = size(params.mask);
  if ~isempty(params.mask) && ((nmask(2)~=nx) || (nmask(1)~=ny))
    if ((nmask(2)==ny) || (nmask(1)==nx))
      error 'Mask array is probably transposed from proper orientation.'
    else
      error 'Mask array must be the same size as the final grid.'
    end
  end
  if ~isempty(params.mask)
    params.maskflag = 1;
  else
    params.maskflag = 0;
  end

  % default for maxiter?
  if isempty(params.maxiter)
    params.maxiter = min(10000,nx*ny);
  end

  % check lengths of the data
  n = length(x);
  if (length(y)~=n) || (length(z)~=n)
    error 'Data vectors are incompatible in size.'
  end
  if n<3
    error 'Insufficient data for surface estimation.'
  end

  % verify the nodes are distinct
  if any(diff(xnodes)<=0) || any(diff(ynodes)<=0)
    error 'xnodes and ynodes must be monotone increasing'
  end
  
  % Are there enough output points to form a surface?
  % Bicubic interpolation requires a 4x4 output grid.  Other types require a 3x3 output grid.
  if strcmp(params.interp, 'bicubic')
    MinAxisLength = 4;
  else
    MinAxisLength = 3;
  end
  if length(xnodes) < MinAxisLength
    error(['The output grid''s x axis must have at least ', num2str(MinAxisLength), ' nodes.']);
  end
  if length(ynodes) < MinAxisLength
    error(['The output grid''s y axis must have at least ', num2str(MinAxisLength), ' nodes.']);
  end
  clear MinAxisLength;

  % do we need to tweak the first or last node in x or y?
  if xmin<xnodes(1)
    switch params.extend
      case 'always'
        xnodes(1) = xmin;
      case 'warning'
%         warning('GRIDFIT:extend',['xnodes(1) was decreased by: ',num2str(xnodes(1)-xmin),', new node = ',num2str(xmin)])
        xnodes(1) = xmin;
      case 'never'
        error(['Some x (',num2str(xmin),') falls below xnodes(1) by: ',num2str(xnodes(1)-xmin)])
    end
  end
  if xmax>xnodes(end)
    switch params.extend
      case 'always'
        xnodes(end) = xmax;
      case 'warning'
%         warning('GRIDFIT:extend',['xnodes(end) was increased by: ',num2str(xmax-xnodes(end)),', new node = ',num2str(xmax)])
        xnodes(end) = xmax;
      case 'never'
        error(['Some x (',num2str(xmax),') falls above xnodes(end) by: ',num2str(xmax-xnodes(end))])
    end
  end
  if ymin<ynodes(1)
    switch params.extend
      case 'always'
        ynodes(1) = ymin;
      case 'warning'
        warning('GRIDFIT:extend',['ynodes(1) was decreased by: ',num2str(ynodes(1)-ymin),', new node = ',num2str(ymin)])
        ynodes(1) = ymin;
      case 'never'
        error(['Some y (',num2str(ymin),') falls below ynodes(1) by: ',num2str(ynodes(1)-ymin)])
    end
  end
  if ymax>ynodes(end)
    switch params.extend
      case 'always'
        ynodes(end) = ymax;
      case 'warning'
        warning('GRIDFIT:extend',['ynodes(end) was increased by: ',num2str(ymax-ynodes(end)),', new node = ',num2str(ymax)])
        ynodes(end) = ymax;
      case 'never'
        error(['Some y (',num2str(ymax),') falls above ynodes(end) by: ',num2str(ymax-ynodes(end))])
    end
  end
  
  % determine which cell in the array each point lies in
  [~, indx] = histc(x,xnodes);
  [~, indy] = histc(y,ynodes);
  % any point falling at the last node is taken to be
  % inside the last cell in x or y.
  k=(indx==nx);
  indx(k)=indx(k)-1;
  k=(indy==ny);
  indy(k)=indy(k)-1;
  ind = indy + ny*(indx-1);
  
  % Do we have a mask to apply?
  if params.maskflag
    % if we do, then we need to ensure that every
    % cell with at least one data point also has at
    % least all of its corners unmasked.
    params.mask(ind) = 1;
    params.mask(ind+1) = 1;
    params.mask(ind+ny) = 1;
    params.mask(ind+ny+1) = 1;
  end
  
  % interpolation equations for each point
  tx = min(1,max(0,(x - xnodes(indx))./dx(indx)));
  ty = min(1,max(0,(y - ynodes(indy))./dy(indy)));
  % Future enhancement: add cubic interpolant
  switch params.interp
    case 'triangle'
      % linear interpolation inside each triangle
      k = (tx > ty);
      L = ones(n,1);
      L(k) = ny;
      
      t1 = min(tx,ty);
      t2 = max(tx,ty);
      A = sparse(repmat((1:n)', 1, 3), [ind, ind + ny + 1, ind + L], [1 - t2, t1, t2 - t1], n, ngrid);
      
    case 'nearest'
      % nearest neighbor interpolation in a cell
      k = round(1-ty) + round(1-tx)*ny;
      A = sparse((1:n)',ind+k,ones(n,1),n,ngrid);
      
    case 'bilinear'
      % bilinear interpolation in a cell
      A = sparse(repmat((1:n)',1,4),[ind,ind+1,ind+ny,ind+ny+1], ...
        [(1-tx).*(1-ty), (1-tx).*ty, tx.*(1-ty), tx.*ty], ...
        n,ngrid);
      
    case 'bicubic'
      % Legacy code calculated the starting index ind for bilinear interpolation, but for bicubic interpolation we need to be further away by one
      % row and one column (but not off the grid).  Bicubic interpolation involves a 4x4 grid of coefficients, and we want x,y to be right
      % in the middle of that 4x4 grid if possible.  Use min and max to ensure we won't exceed matrix dimensions.
      % The sparse matrix format has each column of the sparse matrix A assigned to a unique output grid point.  We need to determine which column
      % numbers are assigned to those 16 grid points.
      % What are the first indexes (in x and y) for the points?
      XIndexes = min(max(1, indx - 1), nx - 3);
      YIndexes = min(max(1, indy - 1), ny - 3);
      % These are the first indexes of that 4x4 grid of nodes where we are doing the interpolation.
      AllColumns = (YIndexes + (XIndexes - 1) * ny)';
      % Add in the next three points.  This gives us output nodes in the first row (i.e. along the x direction).
      AllColumns = [AllColumns; AllColumns + ny; AllColumns + 2 * ny; AllColumns + 3 * ny];
      % Add in the next three rows.  This gives us 16 total output points for each input point.
      AllColumns = [AllColumns; AllColumns + 1; AllColumns + 2; AllColumns + 3];
      % Coefficients are calculated based on:
      % http://en.wikipedia.org/wiki/Lagrange_interpolation
      % Calculate coefficients for this point based on its coordinates as if we were doing cubic interpolation in x.
      % Calculate the first coefficients for x and y.
      XCoefficients = (x(:) - xnodes(XIndexes(:) + 1)) .* (x(:) - xnodes(XIndexes(:) + 2)) .* (x(:) - xnodes(XIndexes(:) + 3)) ./ ((xnodes(XIndexes(:)) - xnodes(XIndexes(:) + 1)) .* (xnodes(XIndexes(:)) - xnodes(XIndexes(:) + 2)) .* (xnodes(XIndexes(:)) - xnodes(XIndexes(:) + 3)));
      YCoefficients = (y(:) - ynodes(YIndexes(:) + 1)) .* (y(:) - ynodes(YIndexes(:) + 2)) .* (y(:) - ynodes(YIndexes(:) + 3)) ./ ((ynodes(YIndexes(:)) - ynodes(YIndexes(:) + 1)) .* (ynodes(YIndexes(:)) - ynodes(YIndexes(:) + 2)) .* (ynodes(YIndexes(:)) - ynodes(YIndexes(:) + 3)));
      % Calculate the second coefficients.
      XCoefficients = [XCoefficients, (x(:) - xnodes(XIndexes(:))) .* (x(:) - xnodes(XIndexes(:) + 2)) .* (x(:) - xnodes(XIndexes(:) + 3)) ./ ((xnodes(XIndexes(:) + 1) - xnodes(XIndexes(:))) .* (xnodes(XIndexes(:) + 1) - xnodes(XIndexes(:) + 2)) .* (xnodes(XIndexes(:) + 1) - xnodes(XIndexes(:) + 3)))];
      YCoefficients = [YCoefficients, (y(:) - ynodes(YIndexes(:))) .* (y(:) - ynodes(YIndexes(:) + 2)) .* (y(:) - ynodes(YIndexes(:) + 3)) ./ ((ynodes(YIndexes(:) + 1) - ynodes(YIndexes(:))) .* (ynodes(YIndexes(:) + 1) - ynodes(YIndexes(:) + 2)) .* (ynodes(YIndexes(:) + 1) - ynodes(YIndexes(:) + 3)))];
      % Calculate the third coefficients.
      XCoefficients = [XCoefficients, (x(:) - xnodes(XIndexes(:))) .* (x(:) - xnodes(XIndexes(:) + 1)) .* (x(:) - xnodes(XIndexes(:) + 3)) ./ ((xnodes(XIndexes(:) + 2) - xnodes(XIndexes(:))) .* (xnodes(XIndexes(:) + 2) - xnodes(XIndexes(:) + 1)) .* (xnodes(XIndexes(:) + 2) - xnodes(XIndexes(:) + 3)))];
      YCoefficients = [YCoefficients, (y(:) - ynodes(YIndexes(:))) .* (y(:) - ynodes(YIndexes(:) + 1)) .* (y(:) - ynodes(YIndexes(:) + 3)) ./ ((ynodes(YIndexes(:) + 2) - ynodes(YIndexes(:))) .* (ynodes(YIndexes(:) + 2) - ynodes(YIndexes(:) + 1)) .* (ynodes(YIndexes(:) + 2) - ynodes(YIndexes(:) + 3)))];
      % Calculate the fourth coefficients.
      XCoefficients = [XCoefficients, (x(:) - xnodes(XIndexes(:))) .* (x(:) - xnodes(XIndexes(:) + 1)) .* (x(:) - xnodes(XIndexes(:) + 2)) ./ ((xnodes(XIndexes(:) + 3) - xnodes(XIndexes(:))) .* (xnodes(XIndexes(:) + 3) - xnodes(XIndexes(:) + 1)) .* (xnodes(XIndexes(:) + 3) - xnodes(XIndexes(:) + 2)))];
      YCoefficients = [YCoefficients, (y(:) - ynodes(YIndexes(:))) .* (y(:) - ynodes(YIndexes(:) + 1)) .* (y(:) - ynodes(YIndexes(:) + 2)) ./ ((ynodes(YIndexes(:) + 3) - ynodes(YIndexes(:))) .* (ynodes(YIndexes(:) + 3) - ynodes(YIndexes(:) + 1)) .* (ynodes(YIndexes(:) + 3) - ynodes(YIndexes(:) + 2)))];
      % Allocate space for all of the data we're about to insert.
      AllCoefficients = zeros(16, n);
      % There may be a clever way to vectorize this, but then the code would be unreadable and difficult to debug or upgrade.
      % The matrix solution process will take far longer than this, so it's not worth the effort to vectorize this.
      for i = 1 : n
        % Multiply the coefficients to accommodate bicubic interpolation.  The resulting matrix is a 4x4 of the interpolation coefficients.
        TheseCoefficients = repmat(XCoefficients(i, :)', 1, 4) .* repmat(YCoefficients(i, :), 4, 1);
        % Add these coefficients to the list.
        AllCoefficients(1 : 16, i) = TheseCoefficients(:);
      end
      % Each input point has 16 interpolation coefficients (because of the 4x4 grid).
      AllRows = repmat(1 : n, 16, 1);
      % Now that we have all of the indexes and coefficients, we can create the sparse matrix of equality conditions.
      A = sparse(AllRows(:), AllColumns(:), AllCoefficients(:), n, ngrid);
  end
  rhs = z;

  % Do we have relative smoothing parameters?
  if numel(params.smoothness) == 1
		% Nothing special; this is just a scalar quantity that needs to be the same for x and y directions.
		xyRelativeStiffness = [1; 1] * params.smoothness;
  else
		% What the user asked for
		xyRelativeStiffness = params.smoothness(:);
  end

  % Build a regularizer using the second derivative.  This used to be called "gradient" even though it uses a second
  % derivative, not a first derivative.  This is an important distinction because "gradient" implies a horizontal
	% surface, which is not correct.  The second derivative favors flatness, especially if you use a large smoothness
	% constant.  Flat and horizontal are two different things, and in this script we are taking an irregular surface and
	% flattening it according to the smoothness constant.
	% The second-derivative calculation is documented here:
	% http://mathformeremortals.wordpress.com/2013/01/12/a-numerical-second-derivative-from-three-points/

	% Minimizes the sum of the squares of the second derivatives (wrt x and y) across the grid
	[i,j] = meshgrid(1:nx,2:(ny-1));
	ind = j(:) + ny*(i(:)-1);
	dy1 = dy(j(:)-1);
	dy2 = dy(j(:));

	Areg = sparse(repmat(ind,1,3),[ind-1,ind,ind+1], ...
		xyRelativeStiffness(2)*[-2./(dy1.*(dy1+dy2)), ...
		2./(dy1.*dy2), -2./(dy2.*(dy1+dy2))],ngrid,ngrid);

	[i,j] = meshgrid(2:(nx-1),1:ny);
	ind = j(:) + ny*(i(:)-1);
	dx1 = dx(i(:) - 1);
	dx2 = dx(i(:));

	Areg = [Areg;sparse(repmat(ind,1,3),[ind-ny,ind,ind+ny], ...
		xyRelativeStiffness(1)*[-2./(dx1.*(dx1+dx2)), ...
		2./(dx1.*dx2), -2./(dx2.*(dx1+dx2))],ngrid,ngrid)];
  nreg = size(Areg, 1);
  
	FidelityEquationCount = size(A, 1);
	% Number of the second derivative equations in the matrix
	RegularizerEquationCount = nx * (ny - 2) + ny * (nx - 2);
	% We are minimizing the sum of squared errors, so adjust the magnitude of the squared errors to make second-derivative
	% squared errors match the fidelity squared errors.  Then multiply by smoothparam.
	NewSmoothnessScale = sqrt(FidelityEquationCount / RegularizerEquationCount);
	
	% Second derivatives scale with z exactly because d^2(K*z) / dx^2 = K * d^2(z) / dx^2.
	% That means we've taken care of the z axis.
	% The square root of the point/derivative ratio takes care of the grid density.
	% We also need to take care of the size of the dataset in x and y.

	% The scaling up to this point applies to local variation.  Local means within a domain of [0, 1] or [10, 11], etc.
	% The smoothing behavior needs to work for datasets that are significantly larger or smaller than that.
	% For example, if x and y span [0 10,000], smoothing local to [0, 1] is insufficient to influence the behavior of
	% the whole surface.  For the same reason there would be a problem applying smoothing for [0, 1] to a small surface
	% spanning [0, 0.01].  Multiplying the smoothing constant by SurfaceDomainScale compensates for this, producing the
	% expected behavior that a smoothing constant of 1 produces noticeable smoothing (when looking at the entire surface
	% profile) and that 1% does not produce noticeable smoothing.
	SurfaceDomainScale = (max(max(xnodes)) - min(min(xnodes))) * (max(max(ynodes)) - min(min(ynodes)));
	NewSmoothnessScale = NewSmoothnessScale *	SurfaceDomainScale;
	
	A = [A; Areg * NewSmoothnessScale];

	rhs = [rhs;zeros(nreg,1)];
  % do we have a mask to apply?
  if params.maskflag
    unmasked = find(params.mask);
  end
  % solve the full system, with regularizer attached
  switch params.solver
    case {'\' 'backslash'}
      if params.maskflag
        % there is a mask to use
        zgrid=nan(ny,nx);
        zgrid(unmasked) = A(:,unmasked)\rhs;
      else
        % no mask
        zgrid = reshape(A\rhs,ny,nx);
      end
      
    case 'normal'
      % The normal equations, solved with \. Can be faster
      % for huge numbers of data points, but reasonably
      % sized grids. The regularizer makes A well conditioned
      % so the normal equations are not a terribly bad thing
      % here.
      if params.maskflag
        % there is a mask to use
        Aunmasked = A(:,unmasked);
        zgrid=nan(ny,nx);
        zgrid(unmasked) = (Aunmasked'*Aunmasked)\(Aunmasked'*rhs);
      else
        zgrid = reshape((A'*A)\(A'*rhs),ny,nx);
      end
      
    case 'symmlq'
      % iterative solver - symmlq - requires a symmetric matrix,
      % so use it to solve the normal equations. No preconditioner.
      tol = abs(max(z)-min(z))*1.e-13;
      if params.maskflag
        % there is a mask to use
        zgrid=nan(ny,nx);
        [zgrid(unmasked),flag] = symmlq(A(:,unmasked)'*A(:,unmasked), ...
          A(:,unmasked)'*rhs,tol,params.maxiter);
      else
        [zgrid,flag] = symmlq(A'*A,A'*rhs,tol,params.maxiter);
        zgrid = reshape(zgrid,ny,nx);
      end
      % display a warning if convergence problems
      switch flag
        case 0
          % no problems with convergence
        case 1
          % SYMMLQ iterated MAXIT times but did not converge.
          warning('GRIDFIT:solver',['Symmlq performed ',num2str(params.maxiter), ...
            ' iterations but did not converge.'])
        case 3
          % SYMMLQ stagnated, successive iterates were the same
          warning('GRIDFIT:solver','Symmlq stagnated without apparent convergence.')
        otherwise
          warning('GRIDFIT:solver',['One of the scalar quantities calculated in',...
            ' symmlq was too small or too large to continue computing.'])
      end
      
    case 'lsqr'
      % iterative solver - lsqr. No preconditioner here.
      tol = abs(max(z)-min(z))*1.e-13;
      if params.maskflag
        % there is a mask to use
        zgrid=nan(ny,nx);
        [zgrid(unmasked),flag] = lsqr(A(:,unmasked),rhs,tol,params.maxiter);
      else
        [zgrid,flag] = lsqr(A,rhs,tol,params.maxiter);
        zgrid = reshape(zgrid,ny,nx);
      end
      
      % display a warning if convergence problems
      switch flag
        case 0
          % no problems with convergence
        case 1
          % lsqr iterated MAXIT times but did not converge.
          warning('GRIDFIT:solver',['Lsqr performed ', ...
            num2str(params.maxiter),' iterations but did not converge.'])
        case 3
          % lsqr stagnated, successive iterates were the same
          warning('GRIDFIT:solver','Lsqr stagnated without apparent convergence.')
        case 4
          warning('GRIDFIT:solver',['One of the scalar quantities calculated in',...
            ' LSQR was too small or too large to continue computing.'])
      end
      
  end  % switch params.solver
  
end  % if params.tilesize...

% only generate xgrid and ygrid if requested.
if nargout>1
  [xgrid,ygrid]=meshgrid(xnodes,ynodes);
end

% ============================================
% End of main function - gridfit
% ============================================

% ============================================
% subfunction - parse_pv_pairs
% ============================================
function params=parse_pv_pairs(params,pv_pairs)
% parse_pv_pairs: parses sets of property value pairs, allows defaults
% usage: params=parse_pv_pairs(default_params,pv_pairs)
%
% arguments: (input)
%  default_params - structure, with one field for every potential
%             property/value pair. Each field will contain the default
%             value for that property. If no default is supplied for a
%             given property, then that field must be empty.
%
%  pv_array - cell array of property/value pairs.
%             Case is ignored when comparing properties to the list
%             of field names. Also, any unambiguous shortening of a
%             field/property name is allowed.
%
% arguments: (output)
%  params   - parameter struct that reflects any updated property/value
%             pairs in the pv_array.
%
% Example usage:
% First, set default values for the parameters. Assume we
% have four parameters that we wish to use optionally in
% the function examplefun.
%
%  - 'viscosity', which will have a default value of 1
%  - 'volume', which will default to 1
%  - 'pie' - which will have default value 3.141592653589793
%  - 'description' - a text field, left empty by default
%
% The first argument to examplefun is one which will always be
% supplied.
%
%   function examplefun(dummyarg1,varargin)
%   params.Viscosity = 1;
%   params.Volume = 1;
%   params.Pie = 3.141592653589793
%
%   params.Description = '';
%   params=parse_pv_pairs(params,varargin);
%   params
%
% Use examplefun, overriding the defaults for 'pie', 'viscosity'
% and 'description'. The 'volume' parameter is left at its default.
%
%   examplefun(rand(10),'vis',10,'pie',3,'Description','Hello world')
%
% params = 
%     Viscosity: 10
%        Volume: 1
%           Pie: 3
%   Description: 'Hello world'
%
% Note that capitalization was ignored, and the property 'viscosity'
% was truncated as supplied. Also note that the order the pairs were
% supplied was arbitrary.

npv = length(pv_pairs);
n = npv/2;

if n~=floor(n)
  error 'Property/value pairs must come in PAIRS.'
end
if n<=0
  % just return the defaults
  return
end

if ~isstruct(params)
  error 'No structure for defaults was supplied'
end

% there was at least one pv pair. process any supplied
propnames = fieldnames(params);
lpropnames = lower(propnames);
for i=1:n
  p_i = lower(pv_pairs{2*i-1});
  v_i = pv_pairs{2*i};
  
  ind = strmatch(p_i,lpropnames,'exact');
  if isempty(ind)
    ind = find(strncmp(p_i,lpropnames,length(p_i)));
    if isempty(ind)
      error(['No matching property found for: ',pv_pairs{2*i-1}])
    elseif length(ind)>1
      error(['Ambiguous property name: ',pv_pairs{2*i-1}])
    end
  end
  p_i = propnames{ind};
  
  % override the corresponding default in params
  params = setfield(params,p_i,v_i); %#ok
  
end


% ============================================
% subfunction - check_params
% ============================================
function params = check_params(params)

% check the parameters for acceptability
% smoothness == 1 by default
if isempty(params.smoothness)
  params.smoothness = 1;
else
  if (numel(params.smoothness)>2) || any(params.smoothness<=0)
    error 'Smoothness must be scalar (or length 2 vector), real, finite, and positive.'
  end
end

% interp must be one of:
% 'bicubic', 'bilinear', 'nearest', or 'triangle'
% but accept any shortening thereof.
valid = {'bicubic', 'bilinear', 'nearest', 'triangle'};
if isempty(params.interp)
  params.interp = 'bilinear';
end
ind = find(strncmpi(params.interp,valid,length(params.interp)));
if (length(ind)==1)
  params.interp = valid{ind};
else
  error(['Invalid interpolation method: ',params.interp])
end

% solver must be one of:
%    'backslash', '\', 'symmlq', 'lsqr', or 'normal'
% but accept any shortening thereof.
valid = {'backslash', '\', 'symmlq', 'lsqr', 'normal'};
if isempty(params.solver)
  params.solver = '\';
end
ind = find(strncmpi(params.solver,valid,length(params.solver)));
if (length(ind)==1)
  params.solver = valid{ind};
else
  error(['Invalid solver option: ',params.solver])
end

% extend must be one of:
%    'never', 'warning', 'always'
% but accept any shortening thereof.
valid = {'never', 'warning', 'always'};
if isempty(params.extend)
  params.extend = 'warning';
end
ind = find(strncmpi(params.extend,valid,length(params.extend)));
if (length(ind)==1)
  params.extend = valid{ind};
else
  error(['Invalid extend option: ',params.extend])
end

% tilesize == inf by default
if isempty(params.tilesize)
  params.tilesize = inf;
elseif (length(params.tilesize)>1) || (params.tilesize<3)
  error 'Tilesize must be scalar and > 0.'
end

% overlap == 0.20 by default
if isempty(params.overlap)
  params.overlap = 0.20;
elseif (length(params.overlap)>1) || (params.overlap<0) || (params.overlap>0.5)
  error 'Overlap must be scalar and 0 < overlap < 1.'
end

% ============================================
% subfunction - tiled_gridfit
% ============================================
function zgrid=tiled_gridfit(x,y,z,xnodes,ynodes,params)
% tiled_gridfit: a tiled version of gridfit, continuous across tile boundaries 
% usage: [zgrid,xgrid,ygrid]=tiled_gridfit(x,y,z,xnodes,ynodes,params)
%
% Tiled_gridfit is used when the total grid is far too large
% to model using a single call to gridfit. While gridfit may take
% only a second or so to build a 100x100 grid, a 2000x2000 grid
% will probably not run at all due to memory problems.
%
% Tiles in the grid with insufficient data (<4 points) will be
% filled with NaNs. Avoid use of too small tiles, especially
% if your data has holes in it that may encompass an entire tile.
%
% A mask may also be applied, in which case tiled_gridfit will
% subdivide the mask into tiles. Note that any boolean mask
% provided is assumed to be the size of the complete grid.
%
% Tiled_gridfit may not be fast on huge grids, but it should run
% as long as you use a reasonable tilesize. 8-)

% Note that we have already verified all parameters in check_params

% Matrix elements in a square tile
tilesize = params.tilesize;
% Size of overlap in terms of matrix elements. Overlaps
% of purely zero cause problems, so force at least two
% elements to overlap.
overlap = max(2,floor(tilesize*params.overlap));

% reset the tilesize for each particular tile to be inf, so
% we will never see a recursive call to tiled_gridfit
Tparams = params;
Tparams.tilesize = inf;

nx = length(xnodes);
ny = length(ynodes);
zgrid = zeros(ny,nx);

% linear ramp for the bilinear interpolation
rampfun = inline('(t-t(1))/(t(end)-t(1))','t');

% loop over each tile in the grid
h = waitbar(0,'Relax and have a cup of JAVA. Its my treat.');
warncount = 0;
xtind = 1:min(nx,tilesize);
while ~isempty(xtind) && (xtind(1)<=nx)
  
  xinterp = ones(1,length(xtind));
  if (xtind(1) ~= 1)
    xinterp(1:overlap) = rampfun(xnodes(xtind(1:overlap)));
  end
  if (xtind(end) ~= nx)
    xinterp((end-overlap+1):end) = 1-rampfun(xnodes(xtind((end-overlap+1):end)));
  end
  
  ytind = 1:min(ny,tilesize);
  while ~isempty(ytind) && (ytind(1)<=ny)
    % update the waitbar
    waitbar((xtind(end)-tilesize)/nx + tilesize*ytind(end)/ny/nx)
    
    yinterp = ones(length(ytind),1);
    if (ytind(1) ~= 1)
      yinterp(1:overlap) = rampfun(ynodes(ytind(1:overlap)));
    end
    if (ytind(end) ~= ny)
      yinterp((end-overlap+1):end) = 1-rampfun(ynodes(ytind((end-overlap+1):end)));
    end
    
    % was a mask supplied?
    if ~isempty(params.mask)
      submask = params.mask(ytind,xtind);
      Tparams.mask = submask;
    end
    
    % extract data that lies in this grid tile
    k = (x>=xnodes(xtind(1))) & (x<=xnodes(xtind(end))) & ...
        (y>=ynodes(ytind(1))) & (y<=ynodes(ytind(end)));
    k = find(k);
    
    if length(k)<4
      if warncount == 0
        warning('GRIDFIT:tiling','A tile was too underpopulated to model. Filled with NaNs.')
      end
      warncount = warncount + 1;
      
      % fill this part of the grid with NaNs
      zgrid(ytind,xtind) = NaN;
      
    else
      % build this tile
      zgtile = RegularizeData3D(x(k),y(k),z(k),xnodes(xtind),ynodes(ytind),Tparams);
      
      % bilinear interpolation (using an outer product)
      interp_coef = yinterp*xinterp;
      
      % accumulate the tile into the complete grid
      zgrid(ytind,xtind) = zgrid(ytind,xtind) + zgtile.*interp_coef;
      
    end
    
    % step to the next tile in y
    if ytind(end)<ny
      ytind = ytind + tilesize - overlap;
      % are we within overlap elements of the edge of the grid?
      if (ytind(end)+max(3,overlap))>=ny
        % extend this tile to the edge
        ytind = ytind(1):ny;
      end
    else
      ytind = ny+1;
    end
    
  end % while loop over y
  
  % step to the next tile in x
  if xtind(end)<nx
    xtind = xtind + tilesize - overlap;
    % are we within overlap elements of the edge of the grid?
    if (xtind(end)+max(3,overlap))>=nx
      % extend this tile to the edge
      xtind = xtind(1):nx;
    end
  else
    xtind = nx+1;
  end

end % while loop over x

% close down the waitbar
close(h)

if warncount>0
  warning('GRIDFIT:tiling',[num2str(warncount),' tiles were underpopulated & filled with NaNs'])
end