function fs=Hypoexponentialpdf(x,l)
fs=0;
for ii=1:length(l)
    prod=1;
    for jj=1:length(l)
        if jj~=ii 
            prod=prod* l(jj)/(l(jj)-l(ii));
        end
    end
    fs=fs+l(ii)*exp(-x*l(ii))*prod;
end
