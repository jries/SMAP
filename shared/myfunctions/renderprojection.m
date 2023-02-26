function img=renderprojection(x,y,z,alpha,beta,rangex,rangey,pixelsize,sigma, ploton,cmax)
%cmax: maximum intensity in rendering: important to keep the intensity
%scaling constant
[xo,yo]=rotcoorddeg(x,y,alpha); % rotation in x,y
[y2,zo]=rotcoorddeg(yo,z,beta); 
rx=rangex(1):pixelsize:rangex(end);
ry=rangey(1):pixelsize:rangey(end);
img=histcounts2(y2,xo,rx,ry);

if sigma>0 %filter
    h=fspecial('gaussian',ceil(2.5*sigma),sigma);
    img=filter2(h,img);
end

if ploton
    figure(88)
    imagesc(img,[0, cmax]);
    colormap hot
    axis equal
end
end


function [xo,yo]=rotcoorddeg(x,y,angle)
xo=cosd(angle)*x+sind(angle)*y;
yo=cosd(angle)*y-sind(angle)*x;
end