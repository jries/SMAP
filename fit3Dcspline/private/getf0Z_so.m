function  f0=getf0Z_so(loc,dz)
%get true bead z-position from the frame in which zfit==0. This is
%determined from a linear fit of zfit vs. frame
%initial position evaluated from z and photons (measure for intensity), as
%sometimes there are erroneous fits with zfit=0/very low photons 
persistent fitt
if isempty(fitt)
    fitt=fittype('poly1');
end
frames=loc.frame;
z=loc.z;
window=round(100/dz);
az=abs(z);
nn=loc.phot./(az+dz);
[~,maxind]=max(nn);

range=max(1,maxind-window):min(maxind+window,length(z));
if sum(abs(z(range))<2*dz)<4
    f0=NaN;
else

    zs=double(z(range));
    fs=double(frames(range));
    fp=fit(zs,fs,fitt);
    f0=fp(0);
end
end