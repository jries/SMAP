function h = fspecial(varargin)
%FSPECIAL Create predefined 2-D filters.
%   H = FSPECIAL(TYPE) creates a two-dimensional filter H of the
%   specified type. Possible values for TYPE are:
%
%     'average'   averaging filter
%     'disk'      circular averaging filter
%     'gaussian'  Gaussian lowpass filter
%     'laplacian' filter approximating the 2-D Laplacian operator
%     'log'       Laplacian of Gaussian filter
%     'motion'    motion filter
%     'prewitt'   Prewitt horizontal edge-emphasizing filter
%     'sobel'     Sobel horizontal edge-emphasizing filter
%     'unsharp'   unsharp contrast enhancement filter
%
%   Depending on TYPE, FSPECIAL may take additional parameters
%   which you can supply.  These parameters all have default
%   values. 
%
%   H = FSPECIAL('average',HSIZE) returns an averaging filter H of size
%   HSIZE. HSIZE can be a vector specifying the number of rows and columns in
%   H or a scalar, in which case H is a square matrix.
%   The default HSIZE is [3 3].
%
%   H = FSPECIAL('disk',RADIUS) returns a circular averaging filter
%   (pillbox) within the square matrix of side 2*RADIUS+1.
%   The default RADIUS is 5.
%
%   H = FSPECIAL('gaussian',HSIZE,SIGMA) returns a rotationally
%   symmetric Gaussian lowpass filter  of size HSIZE with standard
%   deviation SIGMA (positive). HSIZE can be a vector specifying the
%   number of rows and columns in H or a scalar, in which case H is a
%   square matrix.
%   The default HSIZE is [3 3], the default SIGMA is 0.5.
%
%   H = FSPECIAL('laplacian',ALPHA) returns a 3-by-3 filter
%   approximating the shape of the two-dimensional Laplacian
%   operator. The parameter ALPHA controls the shape of the
%   Laplacian and must be in the range 0.0 to 1.0.
%   The default ALPHA is 0.2.
%
%   H = FSPECIAL('log',HSIZE,SIGMA) returns a rotationally symmetric
%   Laplacian of Gaussian filter of size HSIZE with standard deviation
%   SIGMA (positive). HSIZE can be a vector specifying the number of rows
%   and columns in H or a scalar, in which case H is a square matrix.
%   The default HSIZE is [5 5], the default SIGMA is 0.5.
%
%   H = FSPECIAL('motion',LEN,THETA) returns a filter to approximate, once
%   convolved with an image, the linear motion of a camera by LEN pixels,
%   with an angle of THETA degrees in a counter-clockwise direction. The
%   filter becomes a vector for horizontal and vertical motions.  The
%   default LEN is 9, the default THETA is 0, which corresponds to a
%   horizontal motion of 9 pixels.
%
%   H = FSPECIAL('prewitt') returns 3-by-3 filter that emphasizes
%   horizontal edges by approximating a vertical gradient. If you need to
%   emphasize vertical edges, transpose the filter H: H'.
%
%       [1 1 1;0 0 0;-1 -1 -1].
%
%   H = FSPECIAL('sobel') returns 3-by-3 filter that emphasizes
%   horizontal edges utilizing the smoothing effect by approximating a
%   vertical gradient. If you need to emphasize vertical edges, transpose
%   the filter H: H'.
%
%       [1 2 1;0 0 0;-1 -2 -1].
%
%   H = FSPECIAL('unsharp',ALPHA) returns a 3-by-3 unsharp contrast
%   enhancement filter. FSPECIAL creates the unsharp filter from the
%   negative of the Laplacian filter with parameter ALPHA. ALPHA controls
%   the shape of the Laplacian and must be in the range 0.0 to 1.0.
%   The default ALPHA is 0.2.
%
%   Class Support
%   -------------
%   H is of class double.
%
%   Example
%   -------
%      I = imread('cameraman.tif');
%      subplot(2,2,1);imshow(I);title('Original Image'); 
%      H = fspecial('motion',20,45);
%      MotionBlur = imfilter(I,H,'replicate');
%      subplot(2,2,2);imshow(MotionBlur);title('Motion Blurred Image');
%      H = fspecial('disk',10);
%      blurred = imfilter(I,H,'replicate');
%      subplot(2,2,3);imshow(blurred);title('Blurred Image');
%      H = fspecial('unsharp');
%      sharpened = imfilter(I,H,'replicate');
%      subplot(2,2,4);imshow(sharpened);title('Sharpened Image');
%       
%   See also CONV2, EDGE, FILTER2, FSAMP2, FWIND1, FWIND2, IMFILTER.

%   Copyright 1993-2005 The MathWorks, Inc.
%   $Revision: 5.28.4.7 $  $Date: 2006/06/15 20:08:44 $

[type, p2, p3] = ParseInputs(varargin{:});

switch type
  case 'average' % Smoothing filter
     siz = p2;
     h   = ones(siz)/prod(siz);

  case 'disk' % Disk filter
     rad   = p2;
     crad  = ceil(rad-0.5);
     [x,y] = meshgrid(-crad:crad,-crad:crad);
     maxxy = max(abs(x),abs(y));
     minxy = min(abs(x),abs(y));
     m1 = (rad^2 <  (maxxy+0.5).^2 + (minxy-0.5).^2).*(minxy-0.5) + ...
          (rad^2 >= (maxxy+0.5).^2 + (minxy-0.5).^2).* ...
	        sqrt(rad^2 - (maxxy + 0.5).^2);
     m2 = (rad^2 >  (maxxy-0.5).^2 + (minxy+0.5).^2).*(minxy+0.5) + ...
          (rad^2 <= (maxxy-0.5).^2 + (minxy+0.5).^2).* ...
           sqrt(rad^2 - (maxxy - 0.5).^2);
     sgrid = (rad^2*(0.5*(asin(m2/rad) - asin(m1/rad)) + ...
             0.25*(sin(2*asin(m2/rad)) - sin(2*asin(m1/rad)))) - ...
             (maxxy-0.5).*(m2-m1) + (m1-minxy+0.5)) ... 
	          .*((((rad^2 < (maxxy+0.5).^2 + (minxy+0.5).^2) & ...
             (rad^2 > (maxxy-0.5).^2 + (minxy-0.5).^2)) | ...
	          ((minxy==0)&(maxxy-0.5 < rad)&(maxxy+0.5>=rad))));
     sgrid = sgrid + ((maxxy+0.5).^2 + (minxy+0.5).^2 < rad^2);
     sgrid(crad+1,crad+1) = min(pi*rad^2,pi/2);
     if ((crad>0) && (rad > crad-0.5) && (rad^2 < (crad-0.5)^2+0.25)) 
        m1  = sqrt(rad^2 - (crad - 0.5).^2);
	     m1n = m1/rad;
        sg0 = 2*(rad^2*(0.5*asin(m1n) + 0.25*sin(2*asin(m1n)))-m1*(crad-0.5));
        sgrid(2*crad+1,crad+1) = sg0;
        sgrid(crad+1,2*crad+1) = sg0;
        sgrid(crad+1,1)        = sg0;
        sgrid(1,crad+1)        = sg0;
        sgrid(2*crad,crad+1)   = sgrid(2*crad,crad+1) - sg0;
        sgrid(crad+1,2*crad)   = sgrid(crad+1,2*crad) - sg0;
        sgrid(crad+1,2)        = sgrid(crad+1,2)      - sg0;
        sgrid(2,crad+1)        = sgrid(2,crad+1)      - sg0;
     end
     sgrid(crad+1,crad+1) = min(sgrid(crad+1,crad+1),1);
     h = sgrid/sum(sgrid(:));

  case 'gaussian' % Gaussian filter

     siz   = (p2-1)/2;
     std   = p3;
     
     [x,y] = meshgrid(-siz(2):siz(2),-siz(1):siz(1));
     arg   = -(x.*x + y.*y)/(2*std*std);

     h     = exp(arg);
     h(h<eps*max(h(:))) = 0;

     sumh = sum(h(:));
     if sumh ~= 0,
       h  = h/sumh;
     end;
     
  case 'laplacian' % Laplacian filter
     alpha = p2;
     alpha = max(0,min(alpha,1));
     h1    = alpha/(alpha+1); h2 = (1-alpha)/(alpha+1);
     h     = [h1 h2 h1;h2 -4/(alpha+1) h2;h1 h2 h1];

  case 'log' % Laplacian of Gaussian
     % first calculate Gaussian
     siz   = (p2-1)/2;
     std2   = p3^2;
     
     [x,y] = meshgrid(-siz(2):siz(2),-siz(1):siz(1));
     arg   = -(x.*x + y.*y)/(2*std2);
     
     h     = exp(arg);
     h(h<eps*max(h(:))) = 0;

     sumh = sum(h(:));
     if sumh ~= 0,
       h  = h/sumh;
     end;
     % now calculate Laplacian     
     h1 = h.*(x.*x + y.*y - 2*std2)/(std2^2);
     h     = h1 - sum(h1(:))/prod(p2); % make the filter sum to zero
  
  case 'motion' % Motion filter uses bilinear interpolation
     len = max(1,p2);
     half = (len-1)/2;% rotate half length around center
     phi = mod(p3,180)/180*pi;

     cosphi = cos(phi);
     sinphi = sin(phi);
     xsign = sign(cosphi);
     linewdt = 1;

     % define mesh for the half matrix, eps takes care of the right size
     % for 0 & 90 rotation
     sx = fix(half*cosphi + linewdt*xsign - len*eps);
     sy = fix(half*sinphi + linewdt - len*eps);
     [x y] = meshgrid(0:xsign:sx, 0:sy);

     % define shortest distance from a pixel to the rotated line 
     dist2line = (y*cosphi-x*sinphi);% distance perpendicular to the line

     rad = sqrt(x.^2 + y.^2);
     % find points beyond the line's end-point but within the line width
     lastpix = find((rad >= half)&(abs(dist2line)<=linewdt));
     %distance to the line's end-point parallel to the line 
     x2lastpix = half - abs((x(lastpix) + dist2line(lastpix)*sinphi)/cosphi);

     dist2line(lastpix) = sqrt(dist2line(lastpix).^2 + x2lastpix.^2);
     dist2line = linewdt + eps - abs(dist2line);
     dist2line(dist2line<0) = 0;% zero out anything beyond line width

     % unfold half-matrix to the full size
     h = rot90(dist2line,2);
     h(end+(1:end)-1,end+(1:end)-1) = dist2line;
     h = h./(sum(h(:)) + eps*len*len);

     if cosphi>0,
       h = flipud(h);
     end
     
  case 'prewitt' % Prewitt filter
     h = [1 1 1;0 0 0;-1 -1 -1];

  case 'sobel' % Sobel filter
     h = [1 2 1;0 0 0;-1 -2 -1];

  case 'unsharp' % Unsharp filter
     alpha = p2;
     h     = [0 0 0;0 1 0;0 0 0] - fspecial('laplacian',alpha);

  end
  

%%%
%%% ParseInputs
%%%
function [type, p2, p3] = ParseInputs(varargin)

% default values
type      = '';
p2        = [];
p3        = [];

% Check the number of input arguments.
iptchecknargin(1,3,nargin,mfilename);

% Determine filter type from the user supplied string.
type = varargin{1};
type = iptcheckstrs(type,{'gaussian','sobel','prewitt','laplacian','log',...
                    'average','unsharp','disk','motion'},mfilename,'TYPE',1);
  
% default values
switch type
	case 'average'
      p2 = [3 3];  % siz
      
   case 'disk'
      p2 = 5;      % rad
      
   case 'gaussian'
      p2 = [3 3];  % siz
      p3 = 0.5;    % std
      
   case {'laplacian', 'unsharp'}
      p2 = 1/5;    % alpha
      
   case 'log'
      p2 = [5 5];  % siz
      p3 = 0.5;    % std
      
   case 'motion'
      p2 = 9;     % len
      p3 = 0;      % theta
   end
   

switch nargin
    case 1
        % FSPECIAL('average')
        % FSPECIAL('disk')
        % FSPECIAL('gaussian')
        % FSPECIAL('laplacian')
        % FSPECIAL('log')
        % FSPECIAL('motion')
        % FSPECIAL('prewitt')
        % FSPECIAL('sobel')
        % FSPECIAL('unsharp')
        % Nothing to do here; the default values have 
        % already been assigned.        
        
    case 2
       % FSPECIAL('average',N)
       % FSPECIAL('disk',RADIUS)
       % FSPECIAL('gaussian',N)
       % FSPECIAL('laplacian',ALPHA)
       % FSPECIAL('log',N)
       % FSPECIAL('motion',LEN)
       % FSPECIAL('unsharp',ALPHA)
       p2 = varargin{2};
 
       switch type
          case {'sobel','prewitt'}
              msg = sprintf('%s: Too many arguments for this type of filter.', upper(mfilename));
              eid = sprintf('Images:%s:tooManyArgsForThisFilter', mfilename);
              error(eid,msg);
          case {'laplacian','unsharp'}
              iptcheckinput(p2,{'double'},{'nonnegative','real',...
                                  'nonempty','finite','scalar'},...
                            mfilename,'ALPHA',2);
              if  p2 > 1
                  msg = sprintf('%s: ALPHA should be less than or equal 1 and greater than 0.', upper(mfilename));
                  eid = sprintf('Images:%s:outOfRangeAlpha', mfilename);
                  error(eid,msg);
              end
          case {'disk','motion'}
              iptcheckinput(p2,{'double'},...
                            {'positive','finite','real','nonempty','scalar'},...
                            mfilename,'RADIUS or LEN',2);
          case {'gaussian','log','average'}
              iptcheckinput(p2,{'double'},...
                            {'positive','finite','real','nonempty','integer'},...
                            mfilename,'HSIZE',2);
              if numel(p2) > 2
                  msg = 'HSIZE should have 1 or 2 elements.';
                  eid = sprintf('Images:%s:wrongSizeN', mfilename);
                  error(eid,msg);
              elseif numel(p2)==1
                  p2 = [p2 p2]; 
              end
       end       

       
    case 3
       % FSPECIAL('gaussian',N,SIGMA)
       % FSPECIAL('log',N,SIGMA)
       % FSPECIAL('motion',LEN,THETA)
       p2 = varargin{2};
       p3 = varargin{3};
       
       switch type
          case 'motion'
              iptcheckinput(p2,{'double'},...
                            {'positive','finite','real','nonempty','scalar'},...
                            mfilename,'LEN',2);
              iptcheckinput(p3,{'double'},...
                            {'real','nonempty','finite','scalar'},...
                            mfilename,'THETA',3);
          case {'gaussian','log'}
              iptcheckinput(p2,{'double'},...
                            {'positive','finite','real','nonempty','integer'},...
                            mfilename,'N',2);
              iptcheckinput(p3,{'double'},...
                            {'positive','finite','real','nonempty','scalar'},...
                            mfilename,'SIGMA',3);
              if numel(p2) > 2
                  msg = sprintf('%s: size(N) should be less than or equal 2.', upper(mfilename));
                  eid = sprintf('Images:%s:wrongSizeN', mfilename);
                  error(eid,msg);
              elseif numel(p2)==1
                  p2 = [p2 p2]; 
              end
          otherwise   
              msg = sprintf('%s: Too many arguments for this type of filter.', upper(mfilename));
              eid = sprintf('Images:%s:tooManyArgsForThisFilter', mfilename);
              error(eid,msg);
      end
end
