% [o,s]=RWLSPoisson(x,y,N)  : Poisson Reweighted Least Square fit, linear regression with error according to the Poisson distribution. Fits a straight line with offset.
% x : vector of x-values, known postions at which was measured.
% y : vector y-values, measurements.
% N : optional number of measurements from which y was obtained by averaging
% 
% This routine is based on the Wikipedia description at
% https://en.wikipedia.org/wiki/Linear_regression
% (Xt Omega-1 X)^-1  Xt Omega^-1  Y
%
% Example:
% [o,s]=RWLSPoisson([1 2 3 4],[7 8 9 11])
%
function [o,s,vv]=RWLSPoisson(x,y,N)
if nargin < 3
    N=1;
end
myThresh=1.0;  % roughly determined by a simulation

NumIter=5;
v=y; % variances of data is equal (or proportional) to the measured variances
for n=1:NumIter
    vv = v.^2 ./ N; % Variance of the variance. The error of the variance is proportional to the square of the variance, see http://math.stackexchange.com/questions/1015215/standard-error-of-sample-variance
    if any(v<myThresh)
        vv(v<myThresh)=myThresh;  % This is to protect agains ADU-caused bias, which is NOT reduced by averaging
        %if (n==1)
        %    warning('The data has a variance below 2 ADUs at low signal level. This leads to unwanted biases. Increasing the variance estimation for the fit.\n');
        %end
    end
    [o,s]=GLS(x,y,vv);  
    v=o+s*x; % predict variances from the fit for the next round
    %fprintf('RWLSPoisson Iteration %d, o: %g, s: %g\n',n,o,s);
end

