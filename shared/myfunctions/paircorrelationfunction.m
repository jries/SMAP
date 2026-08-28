function [grn,xxo,gr]=paircorrelationfunction(x,y,options)
arguments 
    x
    y
    options.roihandle =[]; %handle to user defined ROI for mask
    options.dr = 2
    options.rmax =[];
end
roi=options.roihandle;
rxrange=[min(x) max(x)];
ryrange=[min(y) max(y)];
edgex=500;
px=options.dr;

rx=rxrange(1)-edgex:px:rxrange(end)+edgex;
if mod(length(rx),2)==1
    rx(end)=[];
end
ry=ryrange(1)-edgex:px:ryrange(end)+edgex;
if mod(length(ry),2)==1
    ry(end)=[];
end

img1=histcounts2(x,y,rx,ry);



%mask
% pm=round(([mean(rx) mean(ry)]-p.sr_pos(1:2))/px);
try
    proi=roi.getPosition;
    mask=roi.createMask;
    % sum(mask,1)f
    roixim=find(sum(mask,1)>0,1,'first');
    roiyim=find(sum(mask,2)>0,1,'first');
    wmaskx=find(sum(mask,1)>0,1,'last')-find(sum(mask,1)>0,1,'first');
    wmasky=find(sum(mask,2)>0,1,'last')-find(sum(mask,2)>0,1,'first');
    switch class(roi)
        case 'imellipse'
            wroi=proi(3:4);
            proimin=proi(1:2);
        otherwise
            wroi =max(proi)-min(proi);
            proimin=min(proi);
    end

    sr_pixrec=(wroi(1)/wmaskx*1000+wroi(2)/wmasky*1000)/2;
    
    maskrs=imresize(mask',sr_pixrec/px);
    dxy=floor((-proimin*1000+[rx(1), ry(1)]+[roixim roiyim]*sr_pixrec)/px);
    mcut=maskrs(dxy(1)+1:dxy(1)+size(img1,1),dxy(2)+1:dxy(2)+size(img1,2));

    % imrgb(:,:,1)=img1;imrgb(:,:,2)=mcut*max(img1(:))/4;imrgb(:,:,3)=0;

    % mp=round(size(maskrs)/2);
    % 
    % rmx=-round(size(img1,1)/2):-round(size(img1,1)/2)+size(img1,1)-1;
    % rmy=-round(size(img1,2)/2):-round(size(img1,2)/2)+size(img1,2)-1;
    % mcut=maskrs(mp(1)+pm(1)+rmx,mp(2)+pm(2)+rmy);
    mcut(mcut>0)=1;
catch err
    mcut=[];
    % err
    % disp('no ROI selected, use FoV instead');
    % mcut=zeros(size(gr1));
    % mcut(1:round(2*p.sr_size(1)/px),1:round(2*p.sr_size(2)/px))=1; %???
    % why? maybe take out edge
end          

[grn,xx,gr]=paircorrelationimage(img1, mcut,options.rmax/options.dr);
xxo=xx*options.dr;
% z=zeros(1001);
% z(1:200:end,:)=1;
% [zgrn,zxx,zgr]=paircorrelationimage(z);

end


function [grn,xx,gr]=paircorrelationimage(img, mask,rmax)
    if nargin<2 || isempty(mask)
        mask=ones(size(img));
    end
    if nargin<3 || isempty(rmax)
        rmax=floor(min(size(img))/2);
    end
    gr1=fftshift(ifft2(  abs(fft2(img)).^2));
    mgr=fftshift(ifft2(  abs(fft2(mask)).^2));
    Amask=sum(mask(:));
    Aimg=sum(img(:));

    rho=Aimg/Amask;

    %normalize
    gr1n=gr1/rho^2./mgr;

    mp=ceil(size(gr1n)/2);

    % nmax=min(mp)-1;
    % rrange=1000;
    % nmax=rrange/px;
    gr1nsm=gr1n(mp(1)-rmax:mp(1)+rmax,mp(2)-rmax:mp(2)+rmax);

    [gr,norm]=radialsum(gr1nsm);
    gr(end)=[];norm(end)=[];
    % nrx=px/2:px:nmax*px+px;
    % xaxnm=nrx;
    % axpc=obj.initaxis('Pair correlation');

    sstart=1;
    
    grn=gr(sstart:end)./norm(sstart:end);
    xx=sstart:length(gr);
    % xaxnms=xaxnm(sstart:end);
end