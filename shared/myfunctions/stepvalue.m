function [sval,nval]=stepvalue(x,istep,fun)
if nargin<3
    fun=@mysimplemean;
end
istep=[istep ;length(x)];
sval=zeros(length(istep)-1,1);
nval=zeros(length(istep)-1,1);
for k=1:length(istep)-1
    sval(k,1)=fun(x(istep(k)+1:istep(k+1)));
    nval(k,1)=istep(k+1)-istep(k);
end
end