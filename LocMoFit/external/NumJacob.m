function df=NumJacob(f,x0,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This code was made by Youngmok Yun, UT Austin. 
% You can distribute or modify as you want, 
% but please do not erase this comment 
% - 2013.05.04
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

epsilon = 1e-6; % delta
l_x0=length(x0); % length of x0;
f0=feval(f,x0',varargin{:}); % caclulate f0
l_f=size(f0,1); % check the size of f
% z = zeros(l_f,1);


for i=1:l_x0
    dx = [ zeros(i-1,1); epsilon; zeros(l_x0-i,1)];
    df(:,i) = ( feval(f,x0'+dx',varargin{:}) - f0)/epsilon;
end
    



