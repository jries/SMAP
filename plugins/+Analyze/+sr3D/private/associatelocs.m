function [beadnum,numlocs]=associatelocs(mx,my,posx,posy,maxd)



beadnum=zeros(length(posx),1);
numlocs=zeros(length(mx),1);
distance2=beadnum+100000000;
for k=1:length(mx)
    r2=(posx-mx(k)).^2+(posy-my(k)).^2;
    ind=r2<maxd^2&r2<distance2;
    distance2(ind)=r2(ind);
    
    beadnum(ind)=k;
    numlocs(k)=sum(ind);
end