function transform=transform_locs_simple(locrefi,loctargeti,p)
df=10;
locref=reducepos(locrefi,df);
loctarget=reducepos(loctargeti,df);

if isfield(p,'tabgroup')
    axh=axes(p.tabgroup);
else
    figure(99);
    axh=gca;
end
axes(axh);
par=axh.Parent;
tg=uitabgroup(par);

tfile=p.Tfile;
 sepscale=2; %maximum separation measure
if exist(tfile,'file')
    l=load(tfile,'transformation');
    Tinitial=l.transformation;
    [loctT.x,loctT.y]=Tinitial.transformCoordinatesInv(loctarget.x(:),loctarget.y(:));
    loctT.frame=loctarget.frame;
    mirrorinfo=Tinitial.tinfo.mirror;
    dx0=0;dy0=0;
%     if contains(mirrorinfo.targetmirror,'no')
%     cutout=true;
%     end
    %     pos=Tinitial.pos;
%     size=Tinitial.size;
else %all initial estimation:
%     approximate shift from size and position
if isfield(p,'separator')
    separator=p.separator;
else
    separator=256;
end
   
    loctT=loctarget;
    locrT=locref;
    separators=[2 2]*separator;
    switch p.Tmode
        case 'up-down'
            loctT.y=loctarget.y(:)-separator;
            targetmirror='no mirror';
            targetpos='bottom';
            separators(2)=separator;
        case 'up-down mirror'
            loctT.y=-loctarget.y(:)+2*separator;
            targetpos='bottom';
            targetmirror='up-down';
            separators(2)=separator;
        case 'right-left'
            loctT.x=loctarget.x(:)-separator;
            targetmirror='no mirror';
            targetpos='right';
            separators(1)=separator;
        case 'right-left mirror'
            loctT.x=-loctarget.x(:)+2*separator;   
            targetmirror='left-right';
            targetpos='right';
            separators(1)=separator;
        case '2 cam'
            %p.Tsplitpos is now the initial scaling factor
            loctT.x=loctarget.x(:)*p.Tsplitpos(1)-min(loctarget.x(:))+10;
            loctT.y=loctarget.y(:)*p.Tsplitpos(end)-min(loctarget.y(:))+10;
            locrT.x=locref.x(:)-min(locref.x(:))+10;
            locrT.y=locref.y(:)-min(locref.y(:))+10;
            separator=max(vertcat(loctT.x(:),loctT.y(:),locrT.x(:),locrT.y(:)))/2;
            sepscale=3;
            targetmirror='no mirror';
            targetpos='center';
    end
    mirrorinfo.targetmirror=targetmirror;
    %determine approximate shift
    xr=1:1:2*separator;yr=xr;
    ht=histcounts2(locrT.x,locrT.y,xr,yr);
    hr=histcounts2(loctT.x,loctT.y,xr,yr);
    G=fftshift(ifft2(conj(fft2(ht)).*fft2(hr)));
    h=fspecial('gaussian',13,2*sepscale);
    Gf=filter2(h,G);
    [~ ,indmax]=max(Gf(:));
    [x0,y0]=ind2sub(size(Gf),indmax);
    dx0=x0-ceil(size(Gf,1)/2);
    dy0=y0-ceil(size(Gf,2)/2);
%     loctT.x=loctT.x-dx;
%     loctT.y=loctT.y-dy;
  

end


[iAa,iBa,na,nb,nseen]=matchlocsall(locrT,loctT,-dx0,-dy0,4*sepscale,1e5);


th=uitab(tg,'Title','intial pos');
axh=axes(th);
plot(axh,locrT.x(iAa),locrT.y(iAa),'ro',loctT.x(iBa)-dx0,loctT.y(iBa)-dy0,'bo',locrT.x(na),locrT.y(na),'r.',loctT.x(nb)-dx0,loctT.y(nb)-dy0,'b.')
transform=interfaces.LocTransform;

t.type=p.Tform;
switch t.type
    case 'polynomial'
        t.parameter=3;
    case 'lwm'
        disp('lwm not implemented')
    case 'pwl'
        disp('pwl not implemented')
end
% t.type='polynomial';
% t.type='projective';
cutofffactor=[1 0.5 0.2 0.1 0.05];
dx=zeros(size(iAa));dy=dx;

% transform.findTransform(locref.x(iAa),locref.y(iAa),loctarget.x(iBa),loctarget.y(iBa),t)
% transform.findTransformZ(locref.x(iAa),locref.y(iAa),locref.z(iAa),loctarget.x(iBa),loctarget.y(iBa),loctarget.z(iBa),t)
% 
%  [xa, ya, za]=transform.transformCoordinatesInv((loctarget.x(iBa)),(loctarget.y(iBa)),(loctarget.z(iBa)));
%     dx=xa-locref.x(iAa);
%    dy=ya-locref.y(iAa);
   
%sort out those which have a large error:
for k=1:length(cutofffactor)
goodind=abs(dx)<cutofffactor(k)*sepscale & abs(dy)<cutofffactor(k)*sepscale;
 transform.findTransform(locref.x(iAa(goodind)),locref.y(iAa(goodind)),loctarget.x(iBa(goodind)),loctarget.y(iBa(goodind)),t)
transform.findTransformZ(locref.x(iAa(goodind)),locref.y(iAa(goodind)),locref.z(iAa(goodind)),loctarget.x(iBa(goodind)),loctarget.y(iBa(goodind)),loctarget.z(iBa(goodind)),t)

 [xa, ya, za]=transform.transformCoordinatesInv((loctarget.x(iBa)),(loctarget.y(iBa)),(loctarget.z(iBa)));
    dx=xa-locref.x(iAa);
   dy=ya-locref.y(iAa);  
end
%    figure(88);plot(locref.x,locref.y,'b.',loctT.x-dx0,loctT.y-dy0,'r+',loctargeti.x,loctargeti.y,'rx',xa,ya,'cx') 
%    

th=uitab(tg,'Title','dx,dy');
axh=axes(th);
dscatter(dx(goodind),dy(goodind))
hold on
circle(0,0,0.02)
axis equal
title(['dx=' num2str(std(dx(goodind)),2) ', dy=' num2str(std(dy(goodind)),2), ', ' num2str(sum(goodind)) ' of ' num2str(nseen) ' paired']);
th=uitab(tg,'Title','beadpos');
axh=axes(th);

 [xaa, yaa, zaa]=transform.transformCoordinatesInv((loctarget.x),(loctarget.y),(loctarget.z));
plot(xaa,yaa,'+',locref.x,locref.y,'o')
legend('target transformed back','reference')

th=uitab(tg,'Title','cross-correlation');
axh=axes(th);
imagesc(Gf);
hold on
plot(y0,x0,'wo',y0,x0,'k+')

transform.tinfo.targetpos=targetpos;
transform.tinfo.separator=separators;
transform.tinfo.mirror=mirrorinfo;
transform.tinfo.cam_pixelsize_nm=1;
transform.tinfo.units='pixels';

end


function pos=reducepos(posin,df)
    z0=ceil(size(posin.x,1)/2);
    for l=size(posin.x,2):-1:1
        framerange=abs(posin.frame(:,l)-z0)<=df;
%         x(:,l)=posin.x(framerange,l);y(:,l)=posin.y(framerange,l);z(:,l)=posin.z(framerange,l);frame(:,l)=posin.frame(framerange,l);
        x(l)=mean(posin.x(framerange,l));y(l)=mean(posin.y(framerange,l));z(l)=mean(posin.z(framerange,l));frame(l)=mean(posin.frame(framerange,l));
    end
    [~,indsort]=sort(frame(:));
%     pos.x=x(indsort);pos.y=y(indsort);pos.z=z(indsort);pos.frame=frame(indsort);
    pos.x=x(indsort)';pos.y=y(indsort)';pos.z=z(indsort)';pos.frame=frame(indsort)';
end