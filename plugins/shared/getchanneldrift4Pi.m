function drift=getchanneldrift4Pi(loc,frameblock,dx0,dy0)
if nargin<3
    dx0=0;dy0=0;
end
frames=[1:frameblock:max(loc.frame) max(loc.frame)+1];
if frames(end)-frames(end-1)<frameblock*0.75 %avoid small last block
    frames(end-1)=[];
end
dx=zeros(length(frames)-1,4);
dy=zeros(length(frames)-1,4);
for k=1:length(frames)-1
    infr=loc.frame>=frames(k) &loc.frame<frames(k+1);
    for ch=2:4
        fnx=['xpix' num2str(ch)];
        dx(k,ch)=getshift(loc.xpix1(infr),loc.(fnx)(infr),loc.xpix1err(infr),loc.([fnx 'err'])(infr));
        fny=['ypix' num2str(ch)];
        dy(k,ch)=getshift(loc.ypix1(infr),loc.(fny)(infr),loc.ypix1err(infr),loc.([fny 'err'])(infr));
    end
end


% for k=1:size(dx,1)
%     dx(k,:)=dx(k,:)-mean(dx(k,:));
%     dy(k,:)=dy(k,:)-mean(dy(k,:));
% end
% x1=mean(dx(:,1));y1=mean(dy(:,1));
% for k=1:size(dx,2) %channel 1: on average shift is zero
%     dx(:,k)=dx(:,k)-x1;
%     dy(:,k)=dy(:,k)-y1;
% end

drift.dxi=dx+dx0;drift.dyi=dy+dy0; drift.framesi=frames;
framesc=(frames(1:end-1)+frames(2:end))/2;
%try to calculate average between channels, so that also ch1 gets drift

%
allframes=(1:max(loc.frame))';
dxs=zeros(length(allframes),4);
dys=zeros(length(allframes),4);
if size(dx,1)<2
    drift.dx=mean(dx,1)+dx0;drift.dy=mean(dy,1)+dy0;drift.frames=1;
    disp(['dx = ' num2str(drift.dx) ', dy = ' num2str(drift.dy)])
else
    for ch=1:4
        dxs(:,ch)=csaps(framesc,dx(:,ch),[],allframes);
        dys(:,ch)=csaps(framesc,dy(:,ch),[],allframes);
    end
    
    figure(334)
    subplot(2,1,1)
    plot(framesc,dx,'o',allframes,dxs)
    subplot(2,1,2)
    plot(framesc,dy,'o',allframes,dys)
    
    drift.dx=dxs+dx0;drift.dy=dys+dy0;drift.frames=allframes;
end
end

function dx=getshift(x1,x2,x1err,x2err)
dxn=x2-x1;
derr2=max(quantile(x1err,0.05).^2+quantile(x2err,0.05).^2,x1err.^2+x2err.^2);
w=1./derr2;
dx=sum(w.*dxn)/sum(w);

% s=median(abs(dxn-dx))/0.6745;
% tune=4.685;
% r=dxn/(tune*s);
% wn=w.*(abs(r)<1) .* (1-r.^2).^2;
% dx=sum(wn.*dxn)/sum(wn);

% goodind=abs(dxn-dx)<2*std(dxn);
% dx=sum(w(goodind).*dxn(goodind))/sum(w(goodind));
% dx=median(dxn);

% dx=mean(dxn);
end
