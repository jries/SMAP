function [srim,nlocs,template]=histrender(pos,rangex, rangey, pixelsx, pixelsy, lut,rangec,template)
% Renders the superresolution image by counting the number of localizations in each pixel
%pos.x, pos.y, pos.c, pos.N,pos.s
%global variables



if isfield(pos,'N')
    N=pos.N;
else
    N=ones(length(pos.x),1,'like',pos.x);
end;

if isfield(pos,'s')
spix=pos.s/pixelsx;
else
    spix=0*pos.x+120;
end
roiks=2.5; %roiks*sigma: size of Roi used in units of sigma
G=creategausstemplate(roiks);
sG=size(G.template);

xpix=(pos.x-rangex(1))/pixelsx;ypix=(pos.y-rangey(1))/pixelsy; %renormalize x, y in units of pixelsize (reconsturcted image)


if nargin<6||length(lut)==1 %only one color
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

[srim,nlocs]=histrenderc(single(xpix-1),single(ypix-1),uint32(srec),single(spix),single(G.template),single(G.sigmatemplate),...
    single(roiks),single(N),int32(uselut),single(pos.c),single(lut), single(rangec));

srim=permute(srim,[2 1 3]);
end

        
function gausstemplate=creategausstemplate(roiks) % create template
% global gausstemplate
% sigmatemplate=10;
sizegauss=5;
sigmatemplate=(sizegauss)/(2*roiks)/2; %for 2.5 sigma in both directions
xg=-sizegauss:sizegauss;
[Xg,Yg]=meshgrid(xg,xg);
template=exp(-((Xg).^2+(Yg).^2)/2/sigmatemplate^2);
gausstemplate.template=template;
gausstemplate.sizegauss=sizegauss;
gausstemplate.sigmatemplate=sigmatemplate;
end
