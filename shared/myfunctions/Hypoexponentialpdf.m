function fs=Hypoexponentialpdf(x,l)
if size(unique(l))<length(l)
l=l+min(l)/100*(0:length(l)-1);
end
fs=0;
for ii=1:length(l)
    prod=1;
    for jj=1:length(l)
        if jj~=ii 
            if l(jj)~=l(ii)
                prod=prod* l(jj)/(l(jj)-l(ii));
            else
                disp('equal k')
            end
               
        end
    end
    fs=fs+l(ii)*exp(-x*l(ii))*prod;
end
