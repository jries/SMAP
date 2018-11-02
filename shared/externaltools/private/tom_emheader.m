function [out] = tom_emheader(in)
%TOM_EMHEADER adds a EM-header structure to a matrix
%
%   [out] = tom_emheader(in)
%
%PARAMETERS
%
%  INPUT
%   in                  (Matrix or Volume)
%  
%  OUTPUT
%   out                 Structure in EM-format with header and in.Value
%   out.Value           Raw data of in
%   out.Header          Header information with standard values
%
%	Build a EM-structure and adds a 512 byte
%   header to the in-Values with default values.
%
%EXAMPLE
%   tom_emheader(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   TOM_EMWRITE, TOM_EMREAD
%
%   created by SN 24.07.2002
%   updated by SN 01.08.2002
%
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2007
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute of Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom

emtype=cellstr(['extern         '; 'EM420          '; 'CM12           '; 'CM200          '; 'CM120/Biofilter'; 'CM300          '; 'Polara         ']);
%                                               EM420=1;CM12=2;CM200=3;
%                                               CM120/Biofilter=4;CM300=5;
%                                               Polara=6;extern=0;

if isequal(computer,'PCWIN') || isequal(computer,'GLNX86') || isequal(computer,'GLNXA64') || isequal(computer,'PCWIN64') || isequal(computer,'MACI')
magic(1)=6; % for PC and Linux 32 and 64 bit
else
magic(1)=3; % for SGI
end;
magic(2)=0;
magic(3)=0;
if isa(in,'double') | isa(in,'single')
magic(4)=5;
end;
if isa(in,'complex') 
magic(4)=8;
end;
if isa(in,'char') 
magic(4)=1;
end;
if isa(in,'int8') 
magic(4)=2;
end;

image_size=size(in)';
comment=char(zeros(80,1));
parameter=zeros(40,1);
fillup=zeros(256,1);

if isstruct(in)~=1
    EM=struct('Magic',magic','Size',image_size,'Comment',comment,'Parameter',parameter,'Fillup',fillup);
Header=struct(...
    'Voltage',parameter(1),...
    'Cs',parameter(2)./1000,...
    'Aperture',parameter(3),...
    'Magnification',parameter(4),...
    'Postmagnification',parameter(5)./1000,...
    'Exposuretime',parameter(6)./1000,...
    'Objectpixelsize',parameter(7)./1000,...
    'Microscope',emtype(parameter(8)+1),...
    'Pixelsize',parameter(9)./1000,...
    'CCDArea',parameter(10)./1000,...
    'Defocus',parameter(11),...
    'Astigmatism',parameter(12),...
    'AstigmatismAngle',parameter(13)./1000,...
    'FocusIncrement',parameter(14)./1000,...
    'CountsPerElectron',parameter(15)./1000,...
    'Intensity',parameter(16)./1000,...
    'EnergySlitwidth',parameter(17),...
    'EnergyOffset',parameter(18),...
    'Tiltangle',parameter(19)./1000,...
    'Tiltaxis',parameter(20)./1000,...
    'Marker_X',parameter(24),...
    'Marker_Y',parameter(25),...
    'Username',num2str(fillup(1:20)),...
    'Date',num2str(fillup(21:28)),...
    'Magic',magic','Size',image_size,'Comment',comment,'Parameter',parameter,'Fillup',fillup,'EM',EM);
    out=struct('Value',in,'Header',Header);

   % Header=struct('Magic',magic','Size',size(in)','Comment',comment,'Parameter',parameter,'Fillup',fillup);
   % out=struct('Value',in,'Header',Header);
else
   disp('Is already an EM-structure !'); 
end;

%
% Structure of EM-Data Files
%
%
%
% Byte 1: Machine Coding:       Machine:    Value:
%                                  OS-9         0
%                                  VAX          1
%                                  Convex       2
%                                  SGI          3
%                                   PC          6
%
% Byte 2: General purpose. On OS-9 system: 0 old version 1 is new version
%
% Byte 3: Not used
%
% Byte 4: Data Type Coding:         Image Type:     No. of Bytes:   Value:
%                                       byte            1               1
%                                       short           2               2
%                                       long int        4               4
%                                       float           4               5
%                                       complex         8               8
%
% Three long integers (3x4 bytes) are image size in x, y, z Dimension
%
% 80 Characters as comment
%
% 40 long integers (4 x 40 bytes) are user defined parameters
%
% Raw data following with the x variable as the fastest dimension, then y and z
%
% 25.6.2002, SN
%



