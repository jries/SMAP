function out=myunique(in)
sin=sort(in);
% out=sin([true; sin(2:end)-sin(1:end-1)>0]);

ind=true(size(sin));
ind(2:end)=diff(sin)>0;

out=sin(ind);

