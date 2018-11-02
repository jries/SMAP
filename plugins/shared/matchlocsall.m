function [iAa,iBa,nA,nB,nseen]=matchlocsall(pos1,pos2,dx,dy,maxd,maxlocsused)
if nargin<6
    maxlocsused=inf;
end


sf1=size(pos1.x);sf2=size(pos2.x);
iAa=zeros(max(sf1(1),sf2(1)),1);
iBa=iAa;

ind1=1; ind2=1;
finda=1;

maxframe=max((pos1.frame(end)),(pos2.frame(end)));

for frame=1:maxframe
    indo1=ind1;indo2=ind2;
    while ind1<=sf1(1)&&(pos1.frame(ind1)<=frame)
        ind1=ind1+1;
    end
    while ind2<=sf2(1)&&(pos2.frame(ind2)<=frame)
        ind2=ind2+1;
    end
    [iA,iB,uiA,uiB]=matchlocs(pos1.x(indo1:ind1-1),pos1.y(indo1:ind1-1),pos2.x(indo2:ind2-1),pos2.y(indo2:ind2-1),([dx dy]),maxd); 
    indi2=finda+length(iA)-1;
    iAa(finda:indi2)=iA+indo1-1;
    iBa(finda:indi2)=iB+indo2-1;
    finda=indi2+1;
    if indi2>maxlocsused
        nseen=ind2;
        break
    end
end
nseen=length(pos2.x);
iAa(finda:end)=[];
iBa(finda:end)=[];

iBtot=1:length(pos2.x);
iAtot=1:length(pos1.x);
nA=setdiff(iAtot,iAa);
nB=setdiff(iBtot,iBa);