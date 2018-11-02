function meanfitpar=meanexp(v,dq,rangev,ax,fac)
q=myquantilefast(abs(v),[0.005,0.995],100000);
if nargin<2||isempty(dq)
    dq=min((q(2)-q(1))/1000,max(max(q(2)/1000,q(1)),q(2)/length(v)*5));
end
if nargin<3 ||isempty(rangev)
    rangev=[min(v) max(v)];
elseif length(rangev)==1
    rangev(2)=max(v);
end

if nargin <5||isempty(fac)
    fac=1;
end
%histogram: center at integer positions
rangev(2)=max(rangev(2),rangev(1)+12*dq);
nh=ceil(rangev(1)):dq:rangev(2);
n=[nh nh(end)+dq]-dq/2;
h=histcounts(v,n);
% nh=n(1:end-1)+(n(2)-n(1))/2;
meanv=num2str(mean(v));
d2=(n(end)-n(1));

fitfun=@(a,b,c,d,e,f,x)(a*exp(b*x)+c*exp(d*x)+e*exp(f*x));
[hx,indx]=max(h);
fo=nh(find(h(indx:end)<hx/2,1,'first')+indx);
fd=7;
f1=3;f2=f1/fd; f3=f2/fd;
startp=[hx*f1^2/3 -f1/fo hx*f2^2/3 -f2/fo hx*f3^2/3 -f3/fo];
cf=fit(nh',h',fitfun,'Lower',[0 -abs(1/n(1)) 0 -abs(1/n(1)) 0 -abs(1/n(1))]*5,...
    'Upper',[Inf -abs(2/n(end)) Inf -abs(2/n(end)) Inf -abs(2/n(end))]*.5,...
    'StartPoint',startp);

% meanfitpar=-(cf.c*cf.b^2+cf.a*cf.d^2)/(cf.c*cf.b^2*cf.d+cf.a*cf.b*cf.d^2);
 meanfitpar2=-(cf.b*cf.d*cf.f*(cf.a/cf.b^2+cf.c/cf.d^2+cf.e/cf.f^2))/(cf.b*cf.d*cf.e+cf.b*cf.c*cf.f+cf.a*cf.d*cf.f);
a=cf.a*exp(cf.b);
c=cf.c*exp(cf.d);
e=cf.e*exp(cf.f);
meanfitpar=(a*(cf.b-1)/cf.b^2+c*(cf.d-1)/cf.d^2+e*(cf.f-1)/cf.f^2)/(a/cf.b+c/cf.d+e/cf.f);
meanfitpar=(meanfitpar+meanfitpar2)/2;

rh=dq/2:dq:rangev(end)+dq;
rh=rangev(1):dq:rangev(end)+dq;
v=cf(rh);
meanfitpar=sum(v.*rh')/sum(v);

if nargin>3&&~isempty(ax)&&ishandle(ax)
%     axis(ax);
%     hold off
% plot(nh,h);
% ax.NextPlot='add';
% hst=fitfun(startp(1),startp(2),startp(3),startp(4),startp(5),startp(6),nh);
% plot(nh,hst,'g')
plot(ax,nh,cf(nh)/fac,'k:');
% hold off
% title(['mfit = ' num2str(meanfitpar) ', mean = ' meanv])
end
