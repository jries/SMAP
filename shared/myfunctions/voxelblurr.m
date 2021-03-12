function v=voxelblurr(fun,par,sigma, pixelsize,rangex, rangey, rangez, varargin)
% fun: function handle
% par: parameters for function
% sigma: scalar, 2-vector (x,y vs z) 
p = inputParser;
p.parameter('fitter',[]);
p.parse(varargin{:})
p = p.Results;
fitter = p.fitter;
if isempty(fitter)
    factor=pixelsize/2; % sampling compared to sigma of Gauss
else
    roiks=2.7; % size of ROI in units of sigma
end
% [x,y,z,norm]=fun(par,single(min(sigma)*factor));
[x,y,z,norm]=fun(par,single(min(sigma)*factor));
indRm = norm == 0;
x = x(~indRm); y = y(~indRm); norm = norm(~indRm);

if ~isempty(z)
    z = z(~indRm);
end
%insert here rotation of coordinates

x=x(:)-rangex(1);y=y(:)-rangey(1);z=z(:)-rangez(1);norm=norm(:);
x=x/pixelsize; y=y/pixelsize; z=z/pixelsize;

srec(1)=ceil((rangex(2)-rangex(1))/pixelsize);
srec(2)=ceil((rangey(2)-rangey(1))/pixelsize);
srec(3)=ceil((rangez(2)-rangez(1))/pixelsize);
sigmax=single(zeros(size(x))+sigma(1));
sigmaz=single(zeros(size(x))+sigma(end));

    if length(z)>0
        [v,nlocs]=gaussrender3c(single(x),single(y),single(z),uint64(srec),...
            single(sigmax),single(sigmaz),...
            single(roiks),single(norm));
    else
        G=creategausstemplate(roiks);
        srec(3) = 1;
        [v,nlocs]=gaussrender_elliptc(single(x),single(y),uint32(srec),single(sigmax),single(sigmax),single(G.template),single(G.sigmatemplate),...
        single(roiks),single(norm),int32(0),single(0),single(0), single([0 0]));
    end
end 

function gausstemplate=creategausstemplate(roiks) % create template
    % global gausstemplate
    % sigmatemplate=10;
    sizegauss=300;
    sigmatemplate=(sizegauss)/(2*roiks)/2; %for 2.5 sigma in both directions
    xg=-sizegauss:sizegauss;
    [Xg,Yg]=meshgrid(xg,xg);
    template=exp(-((Xg).^2+(Yg).^2)/2/sigmatemplate^2);
    gausstemplate.template=template;
    gausstemplate.sizegauss=sizegauss;
    gausstemplate.sigmatemplate=sigmatemplate;
end