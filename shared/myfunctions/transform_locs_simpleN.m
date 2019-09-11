function [transform,iAa,iBa]=transform_locs_simpleN(transform,channelref,locref,channeltarget,loctarget,p)
% locref , target: Mx2 or Mx3 array
if isfield(p,'sepscale')
    sepscale=p.sepscale;
else
sepscale=3;
end



if isfield(p,'Tfile') && exist(p.Tfile,'file')
    Tinitial=load(p.tfile);
    locT=Tinitial.transform2Reference(loctarget);
%     [loctT.x,loctT.y]=Tinitial.transformCoordinatesInv(loctarget.x(:),loctarget.y(:));
%     mirrorinfo=Tinitial.tinfo.mirror;
%     dx=0;dy=0;
else %all initial estimation:
    inforef=transform.info{channelref};
    infotarget=transform.info{channeltarget};
    inforef=takeoutinf(inforef, 2*p.separator);
    infotarget=takeoutinf(infotarget, 2*p.separator);
    locT(:,1)=loctarget(:,1)-infotarget.xrange(1);locT(:,2)=loctarget(:,2)-infotarget.yrange(1);
    locR(:,1)=locref(:,1)-inforef.xrange(1);locR(:,2)=locref(:,2)-inforef.yrange(1);
    
    %mirror if neede
    sref=[inforef.xrange(2)-inforef.xrange(1) inforef.yrange(2)-inforef.yrange(1)];
%     if any(isinf(sref))
%        sref= [2*p.separator 2*p.separator];
%     end
    for k=1:length(inforef.mirror)
        if inforef.mirror(k)>0
            locR(:,inforef.mirror(k))=sref(inforef.mirror(k))-locR(:,inforef.mirror(k));
        end
    end
    star=[infotarget.xrange(2)-infotarget.xrange(1) infotarget.yrange(2)-infotarget.yrange(1)];
    for k=1:length(infotarget.mirror)
        if infotarget.mirror(k)>0
            locT(:,infotarget.mirror(k))=star(infotarget.mirror(k))-locT(:,infotarget.mirror(k));
        end
    end
    

    %determine approximate shift
    xr=1:1:sref(1);yr=1:sref(2);
    ht=histcounts2(locT(:,1),locT(:,2),xr,yr);
    hr=histcounts2(locR(:,1),locR(:,2),xr,yr);
    G=fftshift(ifft2(conj(fft2(hr)).*fft2(ht)));
    h=fspecial('gaussian',13,sepscale);
    Gf=filter2(h,G);
    [~ ,indmax]=max(Gf(:));
    [x0,y0]=ind2sub(size(Gf),indmax);
    dx0=x0-ceil(size(Gf,1)/2);
    dy0=y0-ceil(size(Gf,2)/2);
  

end


locRh.x=locR(:,1);locRh.y=locR(:,2);
locTh.x=locT(:,1);locTh.y=locT(:,2);
if size(locref,2)>3
    locRh.frame=locref(:,4);
    locTh.frame=loctarget(:,4);
else
locRh.frame=ones(size(locRh.x));
locTh.frame=ones(size(locTh.x));
end
[iAa,iBa,na,nb,nseen]=matchlocsall(locRh,locTh,-dx0,-dy0,4*sepscale,1e5);

cutofffactor=[1 0.5 0.2 0.1 0.05 0.03];
dd=zeros(size(iAa,1),2);

for k=1:length(cutofffactor)
goodind=abs(dd(:,1))<cutofffactor(k)*sepscale & abs(dd(:,2))<cutofffactor(k)*sepscale;
 transform.findTransform(channeltarget,locref(iAa(goodind),1:2),loctarget(iBa(goodind),1:2))
 tback=transform.transformToReference(channeltarget,loctarget(iBa,1:2));
dd=tback-locref(iAa,1:2);
end
dd=dd(goodind,:);
% transform.findTransform(channeltarget,locref(iAa,1:2),loctarget(iBa,1:2))


% 
%  tback=transform.transformToReference(channeltarget,loctarget(iBa,1:2));
% dd=tback-locref(iAa,1:2);

%  figure(88);plot(tback(:,1),tback(:,2),'x',locref(:,1),locref(:,2),'o')  
%   figure(88);plot(locref.x,locref.y,'b.',loctT.x,loctT.y,'r+',loctT.x-dx0,loctT.y-dy0,'g.',loctargeti.x,loctargeti.y,'rx',xa,ya,'cx') 
   
if isfield(p,'ax')&& ~isempty(p.ax) && isa(p.ax,'matlab.graphics.axis.Axes')
    axh=p.ax;
    plot(axh,dd(:,1),dd(:,2),'x')
    title(axh,[num2str(std(dd(:,1))) ', ' num2str(std(dd(:,2)))]);
elseif isfield(p,'ax')&& ~isempty(p.ax) && isa(p.ax,'matlab.ui.container.TabGroup') 
     
    axh=axes(uitab(p.ax,'Title','beads'));
    plot(axh,locRh.x,locRh.y,'x',locTh.x,locTh.y,'o')
    plot(axh,locRh.x(iAa),locRh.y(iAa),'ro',locTh.x(iBa)-dx0,locTh.y(iBa)-dy0,'bo')

    
    axh=axes(uitab(p.ax,'Title','difference'));
    plot(axh,dd(:,1),dd(:,2),'x')
    title(axh,[num2str(std(dd(:,1))) ', ' num2str(std(dd(:,2)))]);
    
    axh=axes(uitab(p.ax,'Title','CC'));
    imagesc(axh,Gf);
    axh=axes(uitab(p.ax,'Title','pos'));
    plot(axh,locref(iAa,1),locref(iAa,2),'o',locref(na,1),locref(na,2),'+')
    legend(axh,'paired','not paired')
end
% transform.tinfo.targetpos=targetpos;
% transform.tinfo.separator=separators;
% transform.tinfo.mirror=mirrorinfo;
% transform.tinfo.cam_pixelsize_nm=1;
% transform.tinfo.units='pixels';

end

function struct=takeoutinf(struct, replace)
    struct.xrange(struct.xrange==-Inf)=0;
    struct.xrange(struct.xrange==Inf)=replace;
    struct.yrange(struct.yrange==-Inf)=0;
    struct.yrange(struct.yrange==Inf)=replace;
end

function pos=reducepos(posin,df)
    z0=ceil(size(posin.x,1)/2);
    for l=size(posin.x,2):-1:1
        framerange=abs(posin.frame(:,l)-z0)<=df;
        x(:,l)=posin.x(framerange,l);y(:,l)=posin.y(framerange,l);z(:,l)=posin.z(framerange,l);frame(:,l)=posin.frame(framerange,l);
    end
    [~,indsort]=sort(frame(:));
    pos.x=x(indsort);pos.y=y(indsort);pos.z=z(indsort);pos.frame=frame(indsort);
    
end