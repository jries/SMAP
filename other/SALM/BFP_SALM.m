dna=0.003;
p.n1=1.33;p.n2=1.78;p.lambda=680; p.NA=1.70; p.NAmask=0:dna:1.7;
[r,IS,IU]=intensitySALM(25,p);
dIU=diff(IU);

nax=-p.NA:dna:p.NA;
[X,Y]=meshgrid(nax);
NAXY=sqrt(X.^2+Y.^2);


imout=interp1(p.NAmask(2:end-1)+dna/2,dIU(2:end)./(p.NAmask(2:end-1)+dna/2)/2/pi,NAXY);
% imout=imout/2/pi./NAXY;
figure(88);imagesc(nax,nax,imout)
