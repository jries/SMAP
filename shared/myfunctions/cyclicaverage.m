function midpn=cyclicaverage(dat,period,weights)
if nargin<3
    weights=ones(size(dat));
end
maxerr=mean(dat)*1e-6;

if 0
dat=mod(dat,period);
dat2=dat; 
l=length(dat);
dat2(l+1:2*l)=dat+period;
dat2(2*l+1:3*l)=dat-period;
dat2(3*l+1:4*l)=dat+2*period;
w2=weights;w2(l+1:2*l)=weights;w2(2*l+1:3*l)=weights;w2(3*l+1:4*l)=weights;

% midp=max(dat);
midp=period/2;
%find start
dh=(-0.5:0.2:.5)*period;
for k=length(dh):-1:1
    indh=dat2>midp-period/2+dh(k) & dat2<midp+period/2+dh(k);
%     ssd(k)=std(dat2(indh),w2(indh));  %rms with center?
    ssd(k)=sqrt(mean((dat2(indh)-(midp+dh(k))).^2));
end
% indh=dat2>midp-period/2 & dat2<midp+period/2;
% indh2=dat2>midp & dat2<midp+period;
% indh3=dat2<midp & dat2>midp-period;
% [~, startcase]=max([std(dat2(indh),w2(indh)) std(dat2(indh2),w2(indh2)) std(dat2(indh3),w2(indh3))]);
% startoffs=[-1 0 +1]*period/2;
% [~, startcase]=max(ssd);
[~, startcase]=min(ssd,[],'omitnan');
midp=midp+dh(startcase);

for k=1:50
    indh=dat2>midp-period/2 & dat2<midp+period/2;
    midpn=meanw(dat2(indh),w2(indh));
    if abs(midpn-midp)<maxerr
        break
    end
    midp=midpn;
end
midpn=mod(midpn,period);
end
% alternative approach

ntest=5;
midps=0:period/ntest:period*(1-1/ntest);
for k=length(midps):-1:1
    dath=mod(dat-midps(k)+period/2,period);
    ssdh(k)=sqrt(sum((dath-period/2).^2));
end
[~, startcase]=min(ssdh,[],'omitnan');
midp=midps(startcase);

% midp=0;
error=zeros(50,1);
maxk=20;
for k=1:maxk
    dath=mod(dat-midp+period/2,period);
    midpn=meanw(dath,weights)-period/2+midp;
    error(k)=midpn-midp;
    if abs(midpn-midp)<maxerr
        break
    end
    midp=midpn;
    
end
midpn=mod(midpn,period);
if k==maxk
    disp('cyclic average did not converge')
end





function o=meanw(d,w)
o=sum(d.*w)/sum(w);