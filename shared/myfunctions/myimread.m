function [X, map, alpha] = imread(filename, fmt_s)
%IMREAD Read image from graphics file.
%   A = IMREAD(FILENAME,FMT) reads a grayscale or color image from the file
%   specified by the string FILENAME. If the file is not in the current
%   directory, or in a directory on the MATLAB path, specify the full 
%   pathname.
%   
%   The text string FMT specifies the format of the file by its standard
%   file extension. For example, specify 'gif' for Graphics Interchange 
%   Format files. To see a list of supported formats, with their file 
%   extensions, use the IMFORMATS function. If IMREAD cannot find a file 
%   named FILENAME, it looks for a file named FILENAME.FMT.
%
%   The return value A is an array containing the image data. If the file 
%   contains a grayscale image, A is an M-by-N array. If the file contains
%   a truecolor image, A is an M-by-N-by-3 array. For TIFF files containing
%   color images that use the CMYK color space, A is an M-by-N-by-4 array. 
%   See TIFF in the Format-Specific Information section for more
%   information.
%   
%   The class of A depends on the bits-per-sample of the image data,
%   rounded to the next byte boundary. For example, IMREAD returns 24-bit
%   color data as an array of uint8 data because the sample size for each
%   color component is 8 bits. See the Remarks section for a discussion of  
%   bitdepths, and see the Format-Specific Information section for more  
%   detail about supported bitdepths and sample sizes for a particular
%   format.
%    
%   [X,MAP] = IMREAD(FILENAME,FMT) reads the indexed image in FILENAME into
%   X and its associated colormap into MAP. Colormap values in the image 
%   file are automatically rescaled into the range [0,1]. 
% 
%   [...] = IMREAD(FILENAME) attempts to infer the format of the file
%   from its content.
% 
%   [...] = IMREAD(URL,...) reads the image from an Internet URL.  The
%   URL must include the protocol type (e.g., "http://"). 
%    
%   Remarks
%    
%   Bitdepth is the number of bits used to represent each image pixel.  
%   Bitdepth is calculated by multiplying the bits-per-sample with the 
%   samples-per-pixel. Thus, a format that uses 8-bits for each color 
%   component (or sample) and three samples per pixel has a bitdepth of 24.
%   Sometimes the sample size associated with a bitdepth can be ambiguous: 
%   does a 48-bit bitdepth represent six 8-bit samples or three 16-bit 
%   samples? The following format-specific sections provide sample size 
%   information to avoid this ambiguity.
%    
%   Format-Specific Information (Listed Alphabetically by Format)
%   
%   BMP  --  Windows Bitmap
%
%   Supported  Compression     Output   
%   Bitdepths  None    RLE     Class    Notes
%   ---------------------------------------------------------
%    1-bit      x        -     logical  
%    4-bit      x        x     uint8          
%    8-bit      x        x     uint8
%   16-bit      x        -     uint8    1 sample/pixel
%   24-bit      x        -     uint8    3 samples/pixel
%   32-bit      x        -     uint8    3 samples/pixel (1 byte padding)
%       
%   CUR  -- Cursor File
%  
%   Supported    Compression      Output
%   Bitdepths   None Compressed   Class  
%   --------------------------------------------------
%   1-bit        x      -         logical
%   4-bit        x      -         uint8          
%   8-bit        x      -         uint8
%   
%   Special syntaxes:
%   
%   [...] = IMREAD(...,IDX) reads in one image from a multi-image icon or 
%   cursor file. IDX is an integer value that specifies the order that the
%   image appears in the file. For example, if IDX is 3, IMREAD reads the 
%   third image in the file. If you omit this argument, IMREAD reads the
%   first image in the file. 
% 
%   [A,MAP,ALPHA] = IMREAD(...) returns the AND mask for the resource, 
%   which can be used to determine transparency information.  For cursor 
%   files, this mask may contain the only useful data.    
%     
%   GIF  --  Graphics Interchange Format
%   
%   Supported     Compression        Output
%   Bitdepths    None  Compressed    Class  
%   ---------------------------------------------
%   1-bit         x        -         logical
%   2-to-8 bit    x        -         uint8   
%   
%   Special syntaxes: 
%   
%   [...] = IMREAD(...,IDX) reads in one or more frames from a multiframe 
%   (i.e., animated) GIF file. IDX must be an integer scalar or vector of 
%   integer values.  For example, if IDX is 3, IMREAD reads the third image
%   in the file.  If IDX is 1:5, only the first five frames are returned.
%
%   [...] = IMREAD(...,'Frames',IDX) is the same as the syntax above except
%   that IDX can be 'all'.  In this case, all of the frames are read and 
%   returned in the order that they appear in the file.
%
%   Note: Because of the way GIF files are structured, all of the frames
%   must be read when a particular frame is requested. Consequently, it is 
%   much faster to specify a vector of frames or 'all' for IDX than to call
%   IMREAD in a loop when reading multiple frames from the same GIF file. 
%   
%   HDF  --  Hierarchical Data Format
%     
%   Supported   Raster image   Raster image     Output
%   Bitdepths   with colormap  without colormap Class    Notes
%   ------------------------------------------------------------
%    8-bit        x               x             uint8
%   24-bit        -               x             uint8   3 samples/pixel
%   
%   Special Syntaxes:
%   
%   [...] = IMREAD(...,REF) reads in one image from a multi-image HDF file.
%   REF is an integer value that specifies the reference number used to 
%   identify the image. For example, if REF is 12, IMREAD reads the image 
%   whose reference number is 12. (Note that in an HDF file the reference 
%   numbers do not necessarily correspond with the order of the images in
%   the file. You can use IMFINFO to match up image order with reference 
%   number.) If you omit this argument, IMREAD reads the first image in 
%   the file.
%     
%   ICO  -- Icon File 
%   
%   See CUR.
%   
%   JPEG  --  Joint Photographic Experts Group
%   
%   Note: IMREAD can read any baseline JPEG image as well as JPEG images 
%   with some commonly used extensions. 
%   
%   Supported    Compression      Output
%   Bitdepths   Lossy Lossless    Class      Notes
%   --------------------------------------------------------
%    8-bit        x      x        uint8     Grayscale or RGB
%   12-bit        x      x        uint16    Grayscale  
%   16-bit        -      x        uint16    Grayscale
%   36-bit        x      x        uint16    RGB(Three 12-bit samples/pixel)
%
%   JPEG 2000 - Joint Photographic Experts Group 2000
%
%   Supported      Compression      Output
%   Bitdepths     Lossy Lossless    Class      Notes
%   (per sample)
%   ----------------------------------------------------------
%    1-bit          x      x        logical   Grayscale only
%    2- to 8-bit    x      x        uint8     Grayscale or RGB
%    9- to 16-bit   x      x        uint16    Grayscale or RGB
%
%   Note: Only 1- and 3-sample images are supported.  Indexed JPEG 2000
%   images are not supported.
%
%   Special Syntaxes
%
%   [...] = IMREAD(..., 'Param1', value1, 'Param2', value2, ...) uses
%   parameter-value pairs to control the read operation.  There are two
%   different parameters you can use:
%
%       Parameter name   Value
%       --------------   -----
%       'ReductionLevel' A non-negative integer specifying reduction in  
%                        the resolution of the image. For a reduction 
%                        level 'L', the image resolution is reduced by a 
%                        factor of 2^L. Its default value is 0 implying 
%                        no reduction. The reduction level is limited by 
%                        the total number of decomposition levels as  
%                        provided by 'WaveletDecompositionLevels' field  
%                        in the structure returned from IMFINFO function.   
%
%       'PixelRegion'    {ROWS, COLS}.  IMREAD returns the sub-image
%                        specified by the boundaries in ROWS and COLS.
%                        ROWS and COLS must both be two-element vectors
%                        that denote the 1-based indices [START STOP]. If
%                        'ReductionLevel' is greater than 0, then ROWS and
%                        COLS are coordinates in the reduced-sized image.   
%
%
%   PBM  --  Portable Bitmap
%   
%   Supported  Raw     ASCII (Plain)  Output
%   Bitdepths  Binary  Encoded        Class
%   ----------------------------------------
%   1-bit        x        x          logical
%      
%   PCX  --  Windows Paintbrush
%  
%   Supported     Output    
%   Bitdepths     Class       Notes
%   ----------------------------------------------
%    1-bit        logical     Grayscale only
%    8-bit        uint8       Grayscale or indexed
%   24-bit        uint8       RGB (8-bit samples)
%    
%   PGM  --  Portable Graymap
%        
%   Supported        Raw      ASCII (Plain)  Output        
%   Bitdepths        Binary   Encoded        Class
%   ------------------------------------------------
%   up to 16-bit      x            -         uint8
%   Arbitrary         -            x
%    
%   PNG  --  Portable Network Graphics
%   
%   Supported     Output    
%   Bitdepths     Class      Notes
%   -------------------------------------------
%    1-bit        logical    Grayscale only
%    2-bit        uint8      Grayscale only
%    4-bit        uint8      Grayscale only
%    8-bit        uint8      Grayscale or Indexed
%   16-bit        uint16     Grayscale or Indexed
%   24-bit        uint8      RGB (Three 8-bit samples/pixel)
%   48-bit        uint16     RGB (Three 16-bit samples/pixel)
%         
%   Special Syntaxes:
%   
%   [...] = IMREAD(...,'BackgroundColor',BG) composites any transparent 
%   pixels in the input image against the color specified in BG.  If BG is
%   'none', then no compositing is performed. Otherwise, if the input image
%   is indexed, BG should be an integer in the range [1,P] where P is the
%   colormap length. If the input image is grayscale, BG should be an 
%   integer in the range [0,1].  If the input image is RGB, BG should be a 
%   three-element vector whose values are in the range [0,1]. The string
%   'BackgroundColor' may be abbreviated.  
% 
%   If the ALPHA output argument is used (see below), then BG defaults to 
%   'none' if not specified by the user. Otherwise, if the PNG file 
%   ontains a background color chunk, that color is used as the default  
%   value for BG. If ALPHA is not used and the file does not contain a 
%   background color chunk, then the default value for BG is 1 for indexed  
%   images; 0 for grayscale images; and [0 0 0] for RGB images.  
%
%   [A,MAP,ALPHA] = IMREAD(...) returns the alpha channel if one is
%   present; otherwise ALPHA is [].  Note that MAP may be empty if the file
%   contains a grayscale or truecolor image. 
%     
%   PPM  --  Portable Pixmap 
%   
%   Supported        Raw      ASCII (Plain)  Output        
%   Bitdepths        Binary   Encoded        Class
%   ------------------------------------------------
%   up to 16-bit      x            -         uint8
%   Arbitrary         -            x     
%   
%   RAS  --  Sun Raster 
%   
%   Supported    Output    
%   Bitdepths    Class     Notes
%   ----------------------------------------------------
%    1-bit       logical   Bitmap  
%    8-bit       uint8     Indexed
%   24-bit       uint8     RGB (8-bit samples)
%   32-bit       uint8     RGB with Alpha (8-bit samples)
%    
%   TIFF  --  Tagged Image File Format
%   
%   Supported     Compression                             Output     
%   Bitdepths None Packbits CCITT RGB ICCLAB CIELAB CMYK Class    Notes
%   -----------------------------------------------------------------------
%    1-bit     x      x       x    -    -      -     -   logical  
%    8-bit     x      x       -    -    -      -     -   uint8
%   12-bit     -      -       -    -    -      -     -   uint16  Grayscale
%                                                                or Indexed
%   16-bit     -      -       -    -    -      -     -   uint16  Grayscale 
%                                                                or Indexed
%   24-bit     x      x       -    x    x      x     -   uint8   3 samples
%   32-bit     -      -       -    -    -      -     x   uint8   4 samples
%   36-bit     -      -       -    x    -      -     -   uint16  3 samples
%   48-bit     -      -       -    x    x      x     -   uint16  3 samples      
%   64-bit     -      -       -    -    -      -     x   double  4 samples
%   
%   NOTE: IMREAD supports 8-bit integral and 32-bit floating point tiled 
%   TIFF images, with any compression and colorspace combination listed 
%   above, and 32-bit IEEE floating point images.
%   
%   Special Syntaxes:
%   
%   A = IMREAD(...) returns color data that uses the RGB, CIELAB, ICCLAB,
%   or CMYK color spaces.  If the color image uses the CMYK color space, A 
%   is an M-by-N-by-4 array.
%
%   [...] = IMREAD(..., 'Param1', value1, 'Param2', value2, ...) uses
%   parameter-value pairs to control the read operation.  There are three
%   different parameters you can use:
%
%       Parameter name   Value
%       --------------   -----
%       'Index'          A positive integer specifying which image to read in
%                        a multi-image TIFF file.  For example, if 'Index' is
%                        3, IMREAD reads the third image in the file.
%
%       'Info'           A structure array; the output of IMFINFO.  When
%                        reading images from a multi-image TIFF file, passing
%                        the output of IMFINFO as the 'Info' parameter helps
%                        IMREAD locate the images in the file more quickly.
%
%       'PixelRegion'    {ROWS, COLS}.  IMREAD returns the sub-image
%                        specified by the boundaries in ROWS and COLS.  ROWS
%                        and COLS must be either two- or three-element
%                        vectors.  If two elements are provided, they denote
%                        the 1-based indices [START STOP].  If three elements
%                        are provided, the indices [START INCREMENT STOP]
%                        allow image downsampling.
%   
%   XWD  --  X Window Dump
%   
%   Supported                                  Output    
%   Bitdepths  ZPixmaps  XYBitmaps  XYPixmaps  Class
%   --------------------------------------------------
%   1-bit        x          -         x        logical
%   8-bit        x          -         -        uint8
%
%   Please read the file libtiffcopyright.txt for more information.
%
%   See also IMFINFO, IMWRITE, IMFORMATS, FREAD, IMAGE, DOUBLE, UINT8.

%   Copyright 1984-2008 The MathWorks, Inc.
%   $Revision: 1.1.6.14 $  $Date: 2009/05/14 17:07:43 $


% [filename, fmt_s, extraArgs, msg] = parse_inputs(varargin{:});
% if (~isempty(msg))
%     error('MATLAB:imread:inputParsing', '%s', msg);
% end

% Download remote file.
if (strfind(filename, '://'))
  
    url = true;
  
    if (~usejava('jvm'))
        error('MATLAB:imread:noJava', 'Reading from a URL requires a Java Virtual Machine.')
    end
    
    try
        filename = urlwrite(filename, tempname);
    catch %#ok<*CTCH>
        error('MATLAB:imread:readURL', 'Can''t read URL "%s".', filename);
    end
    
else
  
    url = false;

end

% if (isempty(fmt_s))
%     % The format was not specified explicitly.
%     
%     % Verify that the file exists.
%     fid = fopen(filename, 'r');
%     if (fid == -1)
%       
%         if ~isempty(dir(filename))
%             error('MATLAB:imread:fileOpen', ['Can''t open file "%s" for reading;\nyou' ...
%                    ' may not have read permission.'], ...
%                   filename);
%         else
%             error('MATLAB:imread:fileOpen', 'File "%s" does not exist.', filename);
%         end
% 
%     else
%         % File exists.  Get full filename.
%         filename = fopen(fid);
%         fclose(fid);
%     end
%     
%     % Try to determine the file type.
%     [format, fmt_s] = imftype(filename);
% 
%     if (isempty(format))
%         error('MATLAB:imread:fileFormat', 'Unable to determine the file format.');
%     end
%     
% else
    % The format was specified explicitly.
    
    % Verify that the file exists.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%     fid = fopen(filename, 'r');
%     if (fid == -1)
%         % Couldn't open using the given filename; search for a
%         % file with an appropriate extension.
%         for p = 1:length(fmt_s.ext)
%             fid = fopen([filename '.' fmt_s.ext{p}]);
%             
%             if (fid ~= -1)
%                 % The file was found.  Don't continue searching.
%                 break
%             end
%         end
%     end
%     
%     if (fid == -1)
%         if ~isempty(dir(filename))
%             error('MATLAB:imread:fileOpen', ['Can''t open file "%s" for reading;\nyou' ...
%                    ' may not have read permission.'], ...
%                   filename);
%         else
%             error('MATLAB:imread:fileOpen', 'File "%s" does not exist.', filename);
%         end
%     else
%         filename = fopen(fid);
%         fclose(fid);
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    
% end

if isempty(fmt_s)
    % Get format details.
    fmt_s = imformats(format);
end

% Verify that a read function exists
if (isempty(fmt_s.read))
    error('MATLAB:imread:readFunctionRegistration', 'No reading function for format %s.  See "help imformats".', ...
          fmt_s.ext{1});
end

if ((fmt_s.alpha) && (nargout == 3))
  
    % Use the alpha channel.
    [X, map, alpha] = feval(fmt_s.read, filename, extraArgs{:});
    
else

    % Alpha channel is not requested or is not applicable.
    alpha = [];
%     [X, map] = feval(fmt_s.read, filename, extraArgs{:});
    [X, map] = feval(fmt_s.read, filename);
    
end

% Delete temporary file from Internet download.
if (url)
    delete_download(filename);
end



%%%
%%% Function delete_download
%%%
function delete_download(filename)

try
    delete(filename);
catch
    warning('MATLAB:imread:tempFileDelete', 'Can''t delete temporary file "%s".', filename)
end


      

%%%
%%% Function parse_inputs
%%%
function [filename, fmt_s, extraArgs, msg] = ...
        parse_inputs(varargin)

filename = '';
fmt_s = struct([]);
extraArgs = {};
msg = '';

% Parse arguments based on their number.
switch(nargin)
case 0

    % Not allowed.
    msg = 'Too few input arguments.';
    return;
    
case 1

    % Filename only.
    filename = varargin{1};
    
otherwise

    % Filename and format or other arguments.
    filename = varargin{1};
    
    % Check whether second argument is a format.
    if (ischar(varargin{2}))
        fmt_s = imformats(varargin{2});
    end
    
    if (~isempty(fmt_s))
        % The argument matches a format.
        extraArgs = varargin(3:end);
    else
        % The argument begins the format-specific parameters.
        extraArgs = varargin(2:end);
    end
    
end
