function [img,simulpar]=simulatecamera(locs,p,frames,psf)
%locs.x locs.y locs.z locs.frames locs.phot
%p.xrange p.yrange p.pixelsize p.bg p.psfmode p.calfile
%p.EMon

% p.xrange=[0 20000];p.yrange=[0 20000];p.pixelsize=100;
% p.xrange=[min(locs.x)-1000 max(locs.x)+1000];
% p.yrange=[min(locs.y)-1000 max(locs.y)+1000];
% p.EMon=0;
% p.offset=100; p.conversion=5; p.emgain=100;
% p.background=10;

p.sizex=round((p.xrange(2)-p.xrange(1))/p.pixelsize);
p.sizey=round((p.yrange(2)-p.yrange(1))/p.pixelsize);
im0=zeros(p.sizey,p.sizex,'single');
if isempty(locs)
    simulpar=p;
    img=[];
    return
end
%assume frames sorted, otherwise sort
if nargin<3
    frames=min(locs.frame):max(locs.frame);
end

numlocs=length(locs.x);
%2D in focus
locs.s=zeros(numlocs,1)+100;
img=zeros(p.sizey,p.sizex,length(frames),'single');
ind1=1;

% hwb=waitbar(0,'calcualting camera images');
for k=1:length(frames)
%     waitbar(k/length(frames),hwb);
    while ind1<=numlocs&&locs.frame(ind1)<frames(k)
        ind1=ind1+1;
    end
    if locs.frame(ind1)>frames(k)
        range=[];
    else
        ind2=ind1;
        while ind2<=numlocs&&locs.frame(ind2)==frames(k)
            ind2=ind2+1;
        end    
        range=ind1:ind2-1;
        ind1=ind2;
    end
    
    if ~isempty(range)
        imh=psfimage(locs,range,p,psf);
    else
        imh=im0;
    end
    imh2=addbg(imh,p.background);
    imh3=int2phot(imh2,p.EMon);   
    imh4=phot2adu(imh3,p);
    img(:,:,k)=imh4;
    
    if isfield(p,'plotaxis') &&isvalid(p.plotaxis)
       
        imagesc([ imh4],'Parent',p.plotaxis)
        colorbar('peer',p.plotaxis)
        drawnow
    end
%     waitforbuttonpress;
    
    
end
% close(hwb)
simulpar=p;
end

function imo=phot2adu(imin,p)
if ~p.usecam
    imo=imin;
    return;
end
if ~p.EMon
    emgain=1;
else
    emgain=p.emgain;
end
imo=imin*emgain/p.conversion+p.offset;
end

function img=psfimage(locs,range,p,psf)
locsh.x=locs.x(range);
locsh.y=locs.y(range);
locsh.s=locs.s(range);
locsh.N=locs.phot(range);
locsh.z=locs.znm(range);

% switch p.psfmodel.selection
%     case 'Symmetric Gaussian'
%         [img,nlocs,Gc]=gaussrender(locsh,p.xrange, p.yrange, p.pixelsize, p.pixelsize);        
%     case 'Astigmatig Gaussian'
%     case 'Spline'
        img=psf.render(locsh,p.xrange,p.yrange,p.pixelsize,p.pixelsize);
%         sx=p.sizex;sy=p.sizey;
%         xh=(locsh.x-p.xrange(1))/p.pixelsize;yh=(locsh.y-p.yrange(1))/p.pixelsize;
%         zh=locsh.z;
%         imh=psf.PSF(0,0,0);
%         roipix=size(imh,1);
%         roipixh=floor(roipix/2);
%         
%         
%         imgh=zeros(sx+roipix,sy+roipix);
%         
%         for k=1:length(xh)
%             xr=round(xh(k)); yr=round(yh(k));
%             imh=psf.PSF(xh(k)-xr,yh(k)-yr,zh(k));
%             rxh=xr+1:xr+roipix;
%             ryh=yr+1:yr+roipix;
%             imgh(rxh,ryh)=imgh(rxh,ryh)+imh*locsh.N(k);      
%         end
%         img=imgh(roipixh+1:end-roipixh-1,roipixh+1:end-roipixh-1);
% end

end

function imh2=addbg(imh,bg)
imh2=imh+bg;
end

function imh3=int2phot(imh2,emon)
if emon
    em=2;
else
    em=1;
end
imh3=poissrnd(imh2/em)*em;
end