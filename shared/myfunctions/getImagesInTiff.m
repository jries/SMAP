function ind=getImagesInTiff(name)

file=fopen(name);
bos = fread(file, 2, '*char');
if ( strcmp(bos', 'II') )
    ByteOrder = 'ieee-le';    % Intel little-endian format
elseif ( strcmp(bos','MM') )
    ByteOrder = 'ieee-be';
else
    error('This is not a TIFF file (no MM or II).');
end


%% ---- read in a number which identifies TIFF format

tiff_id = fread(file,1,'uint16', ByteOrder);
ifd_pos   = fread(file, 1, 'uint32', ByteOrder);
stack     = [];
img_indx  = 0;


while  ifd_pos ~= 0 

    img_indx = img_indx + 1;
    fseek(file, ifd_pos, -1);
   

    %read in the number of IFD entries
    num_entries = fread(file,1,'uint16', ByteOrder);
    %fprintf('num_entries = %i\n', num_entries);

    % store current position:
    entry_pos = ifd_pos+2;
    
    % read the next IFD address:
    fseek(file, ifd_pos+12*num_entries+2, -1);
    
    ifd_pos = fread(file, 1, 'uint32', ByteOrder);
end



th=Tiff(name,'r');
th.setDirectory(1);
ind=1;
while ~th.lastDirectory
    ind=ind+1;
    th.setDirectory(ind)
end

% ind=0;
% for lind=16:-1:0
%     if isintiff(th,ind+2^lind)
%         ind=ind+2^lind
%     end
% end

th.close;
% 
% totalstart=2^16;
% start=0;
% max1=getlastintiff(th,start,totalstart)
% maxh=2^ceil(log2(max1)+1);
% 
% for k=1:10
% startnew=getlastintiff(th,start,maxh)
% 
% maxh=2^ceil(log2(startnew-start)+1);
% start=startnew;
% end

% 
% function out=getlastintiff(th,start,max)
% er=true;
% maxh=max;
% while er==true &&maxh>1
%     ind=start+maxh;
%     disp(ind)
%     try
%         th.setDirectory(ind)
%         er=false;
%     catch
%         er=true;
%         maxh=ceil(maxh/2);
%     end
% end
% out=ind;
%     


function out=isintiff(th,num)
try
    th.setDirectory(num)
    out=true;
catch
    out=false;
end