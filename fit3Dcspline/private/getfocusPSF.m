function z0focus=getfocusPSF(PSF)
for k=1:length(PSF)
    [P,CRLB, LL] =mleFit_LM(PSF{k}*10000,4,100,1,0,1);
    sx=P(:,5);
    sy=P(:,6);
    dsx=sx-sy;
    [v,ind(k)]=min(dsx.^2);
end
z0focus=round(mean(ind));
end

