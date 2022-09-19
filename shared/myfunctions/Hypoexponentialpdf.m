function fs=Hypoexponentialpdf(x,l)
fs=0;
for ii=1:length(l)
    prod=1;
    for jj=1:length(l)
        if jj~=ii 
            if l(jj)~=l(ii)
                prod=prod* l(jj)/(l(jj)-l(ii));
            end
               
        end
    end
    fs=fs+l(ii)*exp(-x*l(ii))*prod;
end
