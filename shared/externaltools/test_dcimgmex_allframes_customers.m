%% test mex for dcimg

%% clear all
clc;
clear all; %#ok<CLSCR>
%% open dcimg
% hdcimg        : dcimg handle
% 'FILENAME'    : dcimg file name
FILENAME = 'C:\Hamamatsu\fl4lt_01.dcimg';
hdcimg = dcimgmex( 'open', FILENAME ); %- example
%% get parameter
% param         : return parameter
% hdcimg        : dcimg handle ( returns open )
% 'PARAMNAME'   : parameter name, supported parameters below
%   'SENSOR_BINNING, 'SENSOR_HPOS', 'SENSOR_HSIZE', 'SENSOR_VPOS', 'SENSOR_VSIZE',
%   'IMAGE_WIDTH', 'IMAGE_HEIGHT', 'IMAGE_ROWBYTES', 'IMAGE_PIXELTYPE', 
%   'NUMBEROF_TOTALFRAME', 'NUMBEROF_SESSION', 'NUMBEROF_FRAME', 'NUMBEROF_VIEW', 'NUMBEROF_REGIONRECT',
%   'CURRENT_SESSION", 'CURRENT_VIEW', 'CURRENT_REGIONRECT', 'FILEFORMAT_VERSION'

numFrames = dcimgmex( 'getparam', hdcimg, 'NUMBEROF_FRAME' );

%% read frame
% data          : image data
% hdcimg        : dcimg handle

% Preallocate the array
seq1 = uint16(zeros(2048,2048,1,numFrames)); 
%%seq1(:,:,:,1) = framedatatrans;

for framenum=1:numFrames
   % Read each frame into the appropriate frame in memory.
   data = dcimgmex( 'readframe', hdcimg, framenum);
   tdata = transpose( data );
   seq1(:,:,:,framenum)  = tdata;
   %figure( 1 );  
   %imagesc( tdata );
   %axis off equal;
   colormap gray;
end

montage(seq1);
%implay(seq1,numFrames);

dcimgmex('close', hdcimg);
