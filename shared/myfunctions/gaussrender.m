function [srim,nlocs,Gc]=gaussrender(pos,rangex, rangey, pixelsx, pixelsy, lut,rangec,template)
%pos.x, pos.y, pos.c, pos.N,pos.s
%global variables
persistent G
roiks=2.7;
if isempty(G)
    
    G=creategausstemplate(roiks);
end
Gc=G;

if isfield(pos,'N')
    N=pos.N;
else
    N=ones(length(pos.x),1,'like',pos.x);
end;


% roiks=2.7; %roiks*sigma: size of Roi used in units of sigma
% if nargin==8
%     G=template;
% else
%     G=creategausstemplate(roiks);
% end

% sG=size(G.template);

xpix=(pos.x-rangex(1))/pixelsx;ypix=(pos.y-rangey(1))/pixelsy; %renormalize x, y in units of pixelsize (reconsturcted image)
spix=pos.s/pixelsx;
 spixmed=min(myquantilefast(spix(~isnan(spix)),.5,1000)*10,200);


if nargin<6||length(lut)==1||isempty(lut)||isempty(rangec) %only one color
    uselut=0;
    srec(3)=1;
    lut=0;pos.c=0;rangec=[0,0];
else
    uselut=1;
    srec(3)=3;
%     pos.s(pos.s<1)=.2;
%     pos.c(pos.c>rangec(2))=rangec(2);
%     pos.c(pos.c<rangec(1))=rangec(1);
end





    srec(1)=round((rangex(2)-rangex(1))/pixelsx); %size of reconstructed image. maybe ceil is better?
    srec(2)=round((rangey(2)-rangey(1))/pixelsy);
    
% srim=gaussrenderi(xpix,ypix,srec,pos.s,G.template,G.sigmatemplate,roiks,N,uselut,pos.c,lut, rangec);
if ~isempty(xpix)&&~isempty(ypix)&&~isempty(spix)&&~isempty(pos.c)
    spix(spix>spixmed)=spixmed;
[srim,nlocs]=gaussrenderc(single(xpix-1),single(ypix-1),uint32(srec),single(spix),single(G.template),single(G.sigmatemplate),...
    single(roiks),single(N),int32(uselut),single(pos.c),single(lut), single(rangec));
else
    srim=zeros(srec,'single');
     nlocs=0;
end
srim=permute(srim,[2 1 3]);



end

% 
% figure(3);
% subplot(2,2,2);imagesc(srimc/max(srimc(:)))
% title(sum(srimc(:)))
% colorbar
% subplot(2,2,1);imagesc(srim/max(srim(:)))
% title(max(srim(:)))
% colorbar
% imdiff=srim-srimc;
% subplot(2,2,3);imagesc(imdiff/max(imdiff(:)))
% title(max(imdiff(:)))
% colorbar

% function srim=gaussrenderi(xpix,ypix,srec,sigma,Gtemplate,Gsigma,roiks,N,uselut,c,lut, rangec)
%     s=size(Gtemplate);
%     Gsizegauss=(s(1)-1)/2;
%     sl=length(lut);
%     srim=zeros(srec);
%     for k=1:length(xpix) %all localizations 
%         dn=ceil(roiks*sigma(k));   
%         xr=round(xpix(k));yr=round(ypix(k));
%         dx=xpix(k)-xr;dy=ypix(k)-yr;
%         intcorrection=erf((dn+0.5)/sigma(k)/sqrt(2))^2; %integrate(G,-k sigma, k sigma)= (Erf (k/sqrt(2)))^2: normalization. 0.5: since -dn:dn
%         gaussnorm=N(k)/(2*pi*sigma(k)^2*intcorrection);  
%         if uselut   
%             indc=ceil((c(k)-rangec(1))/(rangec(2)-rangec(1))*(sl));
%         end
% 
%         for xax=-dn:dn
%             xt=round((xax-dx)*Gsigma/sigma(k))+Gsizegauss+1;
%             for yax=-dn:dn
%                 yt=round((yax-dy)*Gsigma/sigma(k))+Gsizegauss+1;
%                 xp=xr+xax;yp=yr+yax;
%                 if xp>0&&xp<=srec(1) && yp>0&&yp<=srec(2)
%                     if uselut
%                         for col=1:3
%                             srim(xp+(yp-1)*srec(1)+(col-1)*srec(1)*srec(2))=srim(xp+(yp-1)*srec(1)+(col-1)*srec(1)*srec(2))...
%                                 +Gtemplate(xt+(yt-1)*sG(1))*gaussnorm*lut(indc+(col-1)*length(lut));
% %                             srim(xp,yp,col)=srim(xp,yp,col)+Gtemplate(xt,yt)*gaussnorm*lut(indc,col);
%                         end
%                     else
%                         srim(xp,yp)=srim(xp,yp)+Gtemplate(xt,yt)*gaussnorm;
%                     end
%                 end
%             end
%         end
% 
%     end
% 
% end
% end

        
function gausstemplate=creategausstemplate(roiks) % create template
% global gausstemplate
% sigmatemplate=10;
sizegauss=600;
sigmatemplate=(sizegauss)/(2*roiks)/2; %for 2.5 sigma in both directions
xg=-sizegauss:sizegauss;
[Xg,Yg]=meshgrid(xg,xg);
template=exp(-((Xg).^2+(Yg).^2)/2/sigmatemplate^2);
gausstemplate.template=template;
gausstemplate.sizegauss=sizegauss;
gausstemplate.sigmatemplate=sigmatemplate;
end
