function ig=withinroi(hroi,x,y,par)

BW = createMask(hroi);

c=par.cpix;
x=round((x-c(3))*par.pixelsize/par.pixrec*1000);
y=round((y-c(1))*par.pixelsize/par.pixrec*1000);
s=size(BW);
x(x>s(1))=s(1);
y(y>s(2))=s(2);
x(x<1)=1;
y(y<1)=1;

lind=sub2ind(size(BW),x,y);

ig=BW(lind);