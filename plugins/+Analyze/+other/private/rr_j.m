function out=rr_j(in)
s=size(in);
nx=(0:s(1)-1)'-s(1)/2;
ny=(0:s(2)-1)'-s(2)/2;
[X,Y]=meshgrid(ny,nx);
out=sqrt(X.^2+Y.^2);