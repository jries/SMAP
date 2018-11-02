function  f0=getf0Z(loc,p)
frames=loc.frame;
z=loc.znm;
window=round(100/p.dz);
az=abs(z);
nn=loc.phot./(az+p.dz);
[~,maxind]=max(nn);

range=max(1,maxind-window):min(maxind+window,length(z));
if sum(abs(z(range))<2*p.dz)<4
    f0=NaN;
else

    zs=double(z(range));
    fs=double(frames(range));
    fp=fit(zs,fs,'poly1');
    f0=fp(0);
end
end