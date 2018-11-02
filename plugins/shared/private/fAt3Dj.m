function pd = fAt3D(zc, yc, xc,xsize,ysize,zsize,delta_f,coeff)
% xc = single(xc);
% yc = single(yc);
% zc = single(zc);
% delta_f = single(delta_f);
% coeff = single(coeff);


% if xc<0
%     xc =0;
% 
% end
% 
% if xc>=xsize
%     xc = xsize-1;
%     
% end
% 
% if yc<0
%     yc =0;
%    
% end
% 
% if yc>=ysize
%     yc = ysize-1;
%     
% end
% 
% if zc<0
%     zc =0;
%     z
% end
% 
% if zc>=zsize
%     zc = zsize-1;
%     
% end
% pd = (0);
%~(xc<0||xc>xsize-1 || yc<0||yc>ysize-1 || zc<0||zc>zsize-1 )

xc = max(xc,0);
xc = min(xc,xsize-1);

yc = max(yc,0);
yc = min(yc,ysize-1);

zc = max(zc,0);
zc = min(zc,zsize-1);


% for i = 0:63
% 
% %     pd = pd+delta_f(i+1)*coeff(xc+1,yc+1,zc+1,i+1);
% 
% %      pd = pd+single(delta_f(i+1)*coeff(i*(xsize*ysize*zsize)+zc*(xsize*ysize)+yc*xsize+xc+1)*xtemp*ytemp*ztemp);
%      pd = pd+(delta_f(i+1)*coeff(i*(xsize*ysize*zsize)+zc*(xsize*ysize)+yc*xsize+xc+1));
% end

pd=sum(delta_f.*(coeff(:,xc+1,yc+1,zc+1)));





