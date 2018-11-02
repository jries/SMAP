function out=gets2z(spline, srange)
% s=p.spline_srange(1): p.spline_ds:p.spline_srange(2);
s=srange;
zall=spline.zrange;
zmat=zeros(length(s));
dmat=zeros(length(s));
h=waitbar(0,'calcualte lookup table');
zold=0;
for k=1:length(s)
    waitbar(k/length(s),h,['calcualte lookup table... ' num2str(k/length(s)*100,'%2.0f') '%']);
    for l=1:length(s)
        [zmat(k,l),dmat(k,l)]=zfromspline(s(k),s(l),spline.x,spline.y,zall,zold);
%         zold=zmat(l,k);
    end
end
delete(h);

z=zall(1):0.01:zall(end);

%interpolation
dfac=5;
ds=(s(2)-s(1))/dfac;
shr=s(1):ds:s(end);
[SX,SY]=meshgrid(s);
[SXhr,SYhr]=meshgrid(shr);
zmathr=interp2(SX,SY,zmat,SXhr,SYhr,'linear');
dmathr=interp2(SX,SY,dmat,SXhr,SYhr,'linear');

% ax=initaxis(p.resultstabgroup,'z');
s2=[shr shr+shr(end)];
dmatp=dmathr/max(dmathr(:))*(max(zmathr(:))-min(zmathr(:)))+min(zmathr(:));
imagesc(s2,s,horzcat(zmathr',dmatp'))
ax=gca;

sx=spline.x(z);sy=spline.y(z);
iin=sx<s(end)&sy<s(end);
hold on;
plot(sx(iin),sy(iin),'w')
plot(sx(iin)+s(end),sy(iin),'w')
hold off;
colorbar(ax)
% ax2=initaxis(p.resultstabgroup,'d');
% imagesc(s,s,dmat)
% colorbar(ax2)
% hold on;
% plot(splinex(z),spliney(z),'w')
% hold off;


% out.z=zmat;
% out.d=dmat;
% out.s=s;
% out.smin=s(1);
% out.smax=s(end);
% out.ds=s(2)-s(1);
% out.zrange=[zall(1) zall(end)];

out.z=single(zmathr);
out.d=single(dmathr);
out.s=single(shr);
out.smin=s(1);
out.smax=s(end);
out.ds=ds;
out.zrange=[zall(1) zall(end)];

[zo,d]=zfromSXSYLut(out,sx,sy);
dz=z(:)-zo(:);
dz(abs(dz)>100)=NaN;
% figure(88);
% plot(z,dz)
out.maxerr=max(abs(dz));
end

function [z,d]=zfromspline(sx,sy,splinex,spliney,zrange,zstart)
% z=fminsearch(@distsplinerr,zstart,[],sx,sy,splinex,spliney,true);
z=fminbnd(@distsplinerr,zrange(1),zrange(end),[],sx,sy,splinex,spliney,true);
 d=distspline(z,sx,sy,splinex,spliney,false);
end

function d=distspline(z,sx,sy,splinex,spliney,sq)
if sq
    d=((sqrt(sx)-sqrt(splinex(z))).^2+(sqrt(sy)-sqrt(spliney(z))).^2);
else
    d=sqrt(((sx)-(splinex(z))).^2+((sy)-(spliney(z))).^2);
end
% 

end
function err= distsplinerr(z,sx,sy,splinex,spliney,sq)
err=sum(distspline(z,sx,sy,splinex,spliney,sq));
end