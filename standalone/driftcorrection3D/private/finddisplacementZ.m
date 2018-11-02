function zpos=finddisplacementZ(xr,zr,xt,zt,xb,zb,window,plotaxis)
if nargin<8
    plotaxis=[];
end
if nargin<7
    window=[];
end


[~,sindr]=sort(xr);
[~,sindt]=sort(xt);
x1r=1;x1t=1;
ccc=zeros(1,length(zb)*2-3);
for k=1:length(xb)-1
   x2r=x1r;x2t=x1t;
   while(xr(sindr(x2r))<xb(k+1))&&x2r<length(sindr)
       x2r=x2r+1;
   end
   while(xt(sindt(x2t))<xb(k+1))&&x2t<length(sindt)
       x2t=x2t+1;
   end
 
   zrh=zr(sindr(x1r:x2r-1));
   zth=zt(sindt(x1t:x2t-1));

   hr=histcounts(zrh,zb);
   ht=histcounts(zth,zb);
   hr=hr-mean(hr);
    ht=ht-mean(ht);
   ch=conv(hr,ht(end:-1:1),'full');
   ccc=ccc+ch;
   x1r=x2r;x1t=x2t; 
end
[mc,ind]=max(ccc);

if isempty(window)
    dh=find(ccc(ind:end)<mc/2,1,'first');
    dh=max(3,round(dh/2));
else
    dh=window;
end
zc=(-length(hr)+1:length(hr)-1)*(zb(2)-zb(1));

inrange=ind-dh:ind+dh;
inrange(zc(inrange)==0)=[];

zred=zc(inrange);
[zpos,fp]=mypeakfit(zc(inrange),ccc(inrange));

indplot=ind-3*dh:ind+3*dh;

if ~isempty(plotaxis) && indplot(1)>0 &&indplot(end)<=length(zc)
  plot(plotaxis,zc(indplot),ccc(indplot),'x')
    plotaxis.NextPlot='add';
    plot(plotaxis,zred,fp(1)*zred.^2+zred*fp(2)+fp(3),'r');
    plotaxis.NextPlot='replace';
    title(plotaxis,['dz (mean): ' num2str(mean(zr)-mean(zt)) ', cc: ' num2str(zpos)])

end


