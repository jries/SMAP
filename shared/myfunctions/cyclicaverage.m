function midpn=cyclicaverage(dat,period,weights)
if nargin<3
    weights=ones(size(dat));
end
dat=mod(dat,period);
maxerr=mean(dat)*1e-6;
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

function o=meanw(d,w)
o=sum(d.*w)/sum(w);