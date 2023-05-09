function hf = NumHessian(f,x0,varargin)
%%variables

epsilon = 1e-5; % delta
l_x0=length(x0); % length of x0;

for i=1:l_x0
    x1 = x0;
    x1(i) = x0(i) - epsilon ;
    df1 = NumJacob(f, x1,varargin{:});
    
    x2 = x0;
    x2(i) = x0(i) + epsilon ;
    df2 = NumJacob(f, x2,varargin{:});
    
    d2f = (df2-df1) / (2*epsilon );
    
    hf(i,:) = d2f';
end
end