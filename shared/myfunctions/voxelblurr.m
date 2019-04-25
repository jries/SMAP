function v=voxelblurr(fun,par,sigma, pixelsize,rangex, rangey, rangez)
% fun: function handle
% par: parameters for function
% sigma: scalar, 2-vector (x,y vs z) 
factor=1; % sampling compared to sigma of Gauss
roiks=2.7; % size of ROI in units of sigma
[x,y,z,norm]=fun(par,single(min(sigma)*factor));

%insert here rotation of coordinates

x=x(:)-rangex(1);y=y(:)-rangey(1);z=z(:)-rangez(1);norm=norm(:);

srec(1)=ceil((rangex(2)-rangex(1))/pixelsize);
srec(2)=ceil((rangey(2)-rangey(1))/pixelsize);
srec(3)=ceil((rangez(2)-rangez(1))/pixelsize);
sigmax=single(zeros(size(x))+sigma(1));
sigmaz=single(zeros(size(x))+sigma(end));

[v,nlocs]=gaussrender3c(single(x),single(y),single(z),uint64(srec),...
    single(sigmax),single(sigmaz),...
    single(roiks),single(norm));

