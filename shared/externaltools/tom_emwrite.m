function tom_emwrite(em_name,Data,form,nr)
% TOM_EMWRITE writes data in an EM-file format
%	 Writes an EM-Image File (V-Format) 
%	 a raw format with a 512 Byte Header.
%	 If input is not a structure a default
%    header is created, which is abandoned 
%    in EMREAD. That way data can be saved 
%    and loaded without providing a header
%    but with compatible file-format to EM.
%
%   tom_emwrite(em_name,Data,form,nr)
%
%PARAMETERS
%
%  INPUT
%   em_name             ['PATHNAME' 'FILENAME']
%   Data                Structure of Image Data
%   Data.Value          Raw data of image, or stack
%   Data.Header         Header information
%  
%  OUTPUT
%
%EXAMPLE
%   tom_emwrite(out);
%      a fileselect-box appears and the data can be saved with the selected
%      filename
%
%      im=tom_emread('pyrodictium_14.em');
%      tom_emwrite('test.em',im);
%
%      modify the comment:
%
%      Data.Header.Comment=double('here is the comment');
%
%REFERENCES
%
%SEE ALSO
%   TOM_EMREAD, TOM_EMHEADER, TOM_READEMHEADER
%
%   created by SN 08/01/02
%   updated by SN 11/07/05 bug fix 'Comment'
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

error(nargchk(1,4,nargin))

  
if nargin <3 form='standard';  nr=1; end;

if nargin <2 form='standard';  nr=1; 
Data=em_name;
[filename, pathname] = uiputfile({'*.em';'*.vol';'*.*'}, 'Save as EM-file');
if isequal(filename,0) | isequal(pathname,0) disp('Data not saved.'); return; end;
em_name=[pathname filename];
if isempty(findstr(em_name,'.'))
    em_name=strcat(em_name,'.em');
end

end;

if nargin <1 error(['Data not specified (e.g. tom_emwrite(out)']);  end;

if isstruct(Data)~=1 & isequal(form,'standard')
    Data=tom_emheader(Data);
end

if isstruct(Data)~=1 & isequal(form,'standard')
    if isequal(computer,'PCWIN') | isequal(computer,'GLNX86') | isequal(computer,'GLNXA64') | isequal(computer,'PCWIN64') | isequal(computer,'MACI')
        magic(1)=6; % for PC and Linux
    else
        magic(1)=3; % for SGI, 
    end;
    magic(2)=0;
    magic(3)=1; % if no structure was given, the user just wants to save in EM-format to be compatible, in emread the header is abandoned. 
    if isa(Data,'double') | isa(Data,'single')
        magic(4)=5;
    end;
    if isa(Data,'complex') 
        magic(4)=8;
    end;
    if isa(Data,'char') 
        magic(4)=1;
    end;
    if isa(Data,'int16') 
        magic(4)=2;
    end;
    comment=char(zeros(80,1));
    parameter=zeros(40,1);
    fillup=zeros(256,1);

    Header=struct('Magic',magic','Size',size(Data)','Comment',comment,'Parameter',parameter,'Fillup',fillup);
    Data=struct('Value',Data,'Header',Header);
    disp('Default EM-header was created.');
end;
emtype=cellstr(['extern         '; 'EM420          '; 'CM12           '; 'CM200          '; 'CM120/Biofilter'; 'CM300          '; 'Polara         '; 'Krios          ']);
%                                               EM420=1;CM12=2;CM200=3;
%                                               CM120/Biofilter=4;CM300=5;
%                                               Polara=6;extern=0;

% set parameter from header structure
Data.Header.EM.Parameter(1)=Data.Header.Voltage;
Data.Header.EM.Parameter(2)=Data.Header.Cs.*1000;
Data.Header.EM.Parameter(3)=Data.Header.Aperture;
Data.Header.EM.Parameter(4)=Data.Header.Magnification;
Data.Header.EM.Parameter(5)=Data.Header.Postmagnification.*1000; 
Data.Header.EM.Parameter(6)=Data.Header.Exposuretime.*1000;
Data.Header.EM.Parameter(7)=Data.Header.Objectpixelsize.*1000;
Data.Header.EM.Parameter(9)=Data.Header.Pixelsize.*1000;
Data.Header.EM.Parameter(10)=Data.Header.CCDArea.*1000;
Data.Header.EM.Parameter(12)=Data.Header.Astigmatism;
Data.Header.EM.Parameter(13)=Data.Header.AstigmatismAngle.*1000;
Data.Header.EM.Parameter(14)=Data.Header.FocusIncrement.*1000;
Data.Header.EM.Parameter(15)=Data.Header.CountsPerElectron.*1000;
Data.Header.EM.Parameter(16)=Data.Header.Intensity.*1000;
Data.Header.EM.Parameter(17)=Data.Header.EnergySlitwidth;
Data.Header.EM.Parameter(18)=Data.Header.EnergyOffset;

Data.Header.EM.Parameter(19)=Data.Header.Tiltangle.*1000;
Data.Header.EM.Parameter(20)=Data.Header.Tiltaxis.*1000;
Data.Header.EM.Parameter(24)=Data.Header.Marker_X;
Data.Header.EM.Parameter(25)=Data.Header.Marker_Y;

%     'Voltage',parameter(1),...
%     'Cs',parameter(2)./1000,...
%     'Aperture',parameter(3),...
%     'Magnification',parameter(4),...
%     'Exposuretime',parameter(6)./1000,...
%     'Objectpixelsize',parameter(7)./1000,...
%     'Microscope',emtype(parameter(8)+1),...
%     'Pixelsize',parameter(9)./1000,...
%     'CCDArea',parameter(10)./1000,...
%     'Defocus',parameter(11),...
%     'Astigmatism',parameter(12),...
%     'AstigmatismAngle',parameter(13)./1000,...
%     'FocusIncrement',parameter(14)./1,...
%     'CountsPerElectron',parameter(15)./1000,...
%     'Intensity',parameter(16)./1000,...
%     'EnergySlitwidth',parameter(17),...
%     'EnergyOffset',parameter(18),...
%     'Tiltangle',parameter(19)./1000,...
%     'Tiltaxis',parameter(20)./1000,...
%     'Username',num2str(fillup(1:20)),...
%     'Date',num2str(fillup(21:28)),...
%     'Magic',magic,'Size',image_size,'Comment',comment,'Parameter',parameter,'Fillup',fillup,'EM',EM);

if isequal(Data.Header.Microscope,'extern') Data.Header.EM.Parameter(8)=0; end;
if isequal(Data.Header.Microscope,'EM420') Data.Header.EM.Parameter(8)=1; end;
if isequal(Data.Header.Microscope,'CM12') Data.Header.EM.Parameter(8)=2; end;
if isequal(Data.Header.Microscope,'CM200') Data.Header.EM.Parameter(8)=3; end;
if isequal(Data.Header.Microscope,'CM120/Biofilter') Data.Header.EM.Parameter(8)=4; end;
if isequal(Data.Header.Microscope,'CM300') Data.Header.EM.Parameter(8)=5; end;
if isequal(Data.Header.Microscope,'Polara') Data.Header.EM.Parameter(8)=6; end;
if isequal(Data.Header.Microscope,'Krios') Data.Header.EM.Parameter(8)=7; end;

Data.Header.EM.Parameter(11)=Data.Header.Defocus;
Data.Header.EM.Parameter(19)=Data.Header.Tiltangle.*1000;
Data.Header.EM.Parameter(20)=Data.Header.Tiltaxis.*1000;
% Comment
%    'Username',num2str(fillup(1:10)),...
%    'Date',num2str(fillup(11:18)),...

% open the stream with the correct format !


    if isequal (computer,'PCWIN') | isequal(computer,'GLNX86') | isequal(computer,'GLNXA64') | isequal(computer,'PCWIN64') | isequal(computer,'MACI')
        fid = fopen(em_name,'w','ieee-le'); Data.Header.Magic(1)=6;
    else
        fid = fopen(em_name,'w','ieee-be');
    end;    
    if fid==-1 error (['Cannot write: ' em_name ' file']); end;
        % writes the header
    %
    % description in 'The Structure of the EM-Data Files', Herr Hegerl
    %

    
    
    if size(Data.Header.Magic,1)~=4
    	Data.Header.Magic=[6 0 0 5]';
    end
    fwrite(fid,Data.Header.Magic,'char');

    if size(Data.Header.Size,1)~=3
	Data.Header.Size(1)=size(Data.Value,1);
	Data.Header.Size(2)=size(Data.Value,2);
	Data.Header.Size(3)=size(Data.Value,3);
    end
    fwrite(fid,Data.Header.Size,'int32');
    if size(Data.Header.Comment,2)<80
        f=ones((80-size(Data.Header.Comment,2)),1).*32;
    	Data.Header.Comment(size(Data.Header.Comment,2)+1:80)=f;
    else
        Data.Header.Comment=Data.Header.Comment(1:80);
    end
    fwrite(fid,Data.Header.Comment,'char');

    if size(Data.Header.EM.Parameter,1)~=40
    	Data.Header.EM.Parameter((end+1):40)=ones(40-length(Data.Header.EM.Parameter),1);
    end
    fwrite(fid,Data.Header.EM.Parameter,'int32');

    if (~isfield(Data.Header,'Fillup') || size(Data.Header.Fillup,1)~=256) 
    	Data.Header.Fillup=ones(256,1);
    end
fwrite(fid,Data.Header.Fillup,'char');

% the size of the image
% to adapt to EM, transpose '
%xdim = Data.Header.Size(2);
%ydim = Data.Header.Size(1);

xdim = Data.Header.Size(1);
ydim = Data.Header.Size(2);
zdim = Data.Header.Size(3);
Data_write=0;

for lauf=1:zdim
	% what byte format?
	% to adapt to EM, transpose '
	%Data_write=Data.Value(1:xdim,1:ydim,lauf)';
	
	Data_write=squeeze(Data.Value(1:xdim,1:ydim,lauf));

	if Data.Header.Magic(4)==1
		fwrite(fid,Data_write,'char');
	elseif Data.Header.Magic(4)==2
		fwrite(fid,Data_write,'int16');
	elseif Data.Header.Magic(4)==4
		fwrite(fid,Data_write,'int32');
	elseif Data.Header.Magic(4)==5
		fwrite(fid,Data_write,'float');
	elseif Data.Header.Magic(4)==8
		fwrite(fid,Data_write,'float64');
	else
		disp('Sorry, i cannot write this as an EM-File !!!');
	end;

end;

fclose(fid);

