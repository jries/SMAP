function fs=Erlangdistribution(x,l,k)
fs=l^k*x.^(k-1).*exp(-l*x)/factorial(k-1);
end