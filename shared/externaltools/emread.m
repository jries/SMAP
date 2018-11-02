function I = emread(em_name)
% emread -- EM to Matlab cube converter
%
%  Usage:
%    I = emread(em_name)
%  Inputs:
%    em_name    3-d cube in EM format (V, floating point / byte)
%
%  Outputs:
%    I    	cube
%
%  See Also
%    emwrite, emreadheader
%

%    Last changes:
%    Oct. 27, 2003
%    M. Riedlberger


%f=sprintf('Reading EM-file: %s',em_name);disp(f);

fid = fopen(em_name,'r','ieee-le');
if fid==-1 error(sprintf('Wrong File Name: %s',em_name)); end;
machine = fread(fid,[1],'uint8');
header = fread(fid,[2],'uint8');
data_type = fread(fid,[1],'uint8');
fclose(fid);
%disp(machine)
if machine==6
    fid = fopen(em_name,'r','ieee-le');
elseif machine==3
    fid = fopen(em_name,'r','ieee-be');
elseif machine==5
    fid = fopen(em_name,'r','ieee-be');
else
    error('Error: Wrong File Format');
end

header = fread(fid,[128],'uint32');
xdim = header(2);
ydim = header(3);
zdim = header(4);
%f=sprintf('Reading EM-file: %s with Dimensions:x=%g, y=%g, z=%g',em_name, xdim,ydim,zdim); disp(f);

%disp(data_type)

if data_type==5
    I = reshape(fread(fid,[xdim * ydim * zdim],'float'),xdim, ydim, zdim);
elseif data_type==1
    I = uint8(reshape(fread(fid,[xdim * ydim * zdim],'uint8'),xdim, ydim, zdim));
elseif data_type==2
    I = uint8(reshape(fread(fid,[xdim * ydim * zdim],'uint16'),xdim, ydim, zdim));
else
    error('Error: Wrong Data Type');
end

fclose(fid);
