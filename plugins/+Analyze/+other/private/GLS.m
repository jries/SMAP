% [o,s]=GLS(x,y,v)  : Generalized Least Square fit, linear regression with error. Fits a straight line with offset through data with known variances.
% x : vector of x-values, known postions at which was measured.
% y : vector y-values, measurements.
% v : vector of known errors. These can also be estimated, but if you have a model, use RWLS, the reweighted least squares regression.
%     Specifically look at RWLSPoisson if this fitting is for Poisson statistics.
% o : offset = y-value at zero x-value
% s : slope
% 
% This routine is based on the Wikipedia description at
% https://en.wikipedia.org/wiki/Linear_regression
% (Xt Omega-1 X)^-1  Xt Omega^-1  Y
%
% Example:
% [o,s]=GLS([1 2 3 4],[7 8 9 11],[1 1 1 1])
%
function [o,s]=GLS(x,y,v)

X=transpose([ones(1,size(x,2));x]);
Xt=transpose(X);
Omega_X=transpose([1./v;x./v]);
Omega_Y=transpose(y./v);
ToInvert = (Xt * Omega_X);

res=inv(ToInvert)*(Xt*Omega_Y); 

o=res(1);
s=res(2);
