%make rotated Gaussian PSF

x1 = -3000:100:3000; x2 = -3000:100:3000;
[X1,X2] = meshgrid(x1, x2);
sigmax0=100;
sigmay0=130;
d=250;
lambda=600;
n=1.33;
Theta = 15;

zRx=pi*sigmax0^2/lambda*n*4;
zRy=pi*sigmay0^2/lambda*n*4;

z=-1000:10:1000;



outim=zeros(length(x1),length(x2),length(z));
for k=1:length(z)
sigma1=sigmax0*sqrt(1+((z(k)-d)/zRx).^2);
sigma2=sigmay0*sqrt(1+((z(k)+d)/zRy).^2);
a = ((cosd(Theta)^2) / (2*sigma1^2)) + ((sind(Theta)^2) / (2*sigma2^2));
b = -((sind(2*Theta)) / (4*sigma1^2)) + ((sind(2*Theta)) / (4*sigma2^2));
c = ((sind(Theta)^2) / (2*sigma1^2)) + ((cosd(Theta)^2) / (2*sigma2^2));

mu = [0 0];
A = pi/sigma1/sigma2;
outim(:,:,k) = A*exp(-(a*(X1 - mu(1)).^2 + 2*b*(X1 - mu(1)).*(X2 - mu(2)) + c*(X2 - mu(2)).^2));
end

imx(outim)