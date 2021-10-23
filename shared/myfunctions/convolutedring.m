function out=convolutedring(p,xv,sigma)
%p=[A,x0,R], 
A=p(1);
x0=p(2);
R=p(3);
off=p(4);
out=quadv(@(x)convring(x,xv,sigma,R,x0),-R+x0,R+x0);
out=A*out/max(out)+off;
end

function Y = convring(xs,x,pp,ps,x0)
Y=profileh(x-xs,pp).*structurering(xs-x0,ps);
end

function out=profileh(x,pp)
out=exp(-x.^2/2/pp(1)^2);
end

function out=structurering(x,ps)
out=ps(1)/sqrt(ps(1)^2-x.^2);
if imag(out)
    out=-1000+0*x;
end
end