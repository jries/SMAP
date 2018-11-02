function out=runningaverage(in,dn,varargin)
if ~iscell(in) %to work also with vectors
    in=num2cell(in);
    convert=true;
else
    convert=false;
end
out{length(in)}=0*in{1};
for k=1:length(in)
    out{k}=0*in{k};
    l1=round(max(1,k-(dn-1)/2));
    l2=round(min(length(in),k+(dn-1)/2));
    for l=l1:l2
         im1=out{k};
        im2=in{l};
%          whos im1 im2
        out{k}=im1+im2;
    end
    norm=1/(l2-l1+1);
    out{k}=out{k}*norm;
end
if convert
    out=cell2mat(out);
end
