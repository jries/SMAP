%Simulate kinetic processes, chain of events with kinetic constants ki

N=1000000; %number of data points
tausingle=[10 30];
eps=0.01; % ki ~= kj for all i,j
taus=[tausingle, tausingle+eps];
% taus=[10.001 30.0 10.0 30.001]; % list of kinetic constants
dt=1;

timet=zeros(N,length(ks));

for k=1:length(taus)
    timet(:,k)=exprnd(taus(k),N,1);
end
ks=1./taus;

%total cycle time:
t2=sum(timet,2);
figure(88);
subplot(2,1,1)
[fp,nt,ht]=plotwithfit(t2,dt);
title(['16 nm steps, tau=' num2str(1/fp.k1)])
fx=Hypoexponentialpdf(nt,ks);
plot(nt,fx*N*dt,'r')

%try fit
% fitf=fittype('N*dt*Hypoexp4fit(x,k1,k2)');
fitf=@(k1,k2,A,x) N*dt*Hypoexp4fit(x,k1,k2,A);

[htmax,imax]=max(ht);
kstart=1/nt(imax);
plot(nt,fitf(kstart*2,kstart*4,1,nt),'g--')
fph=fit(nt',ht',fitf,'StartPoint',[kstart*2,kstart*4,1]);
plot(nt,fph(nt),'b')


subplot(2,1,2)
l=size(timet,2); 
midp=round(l/2);
t8_1=sum(timet(:,1:midp),2);
t8_2=sum(timet(:,end-midp+1:end),2);
fp=plotwithfit(vertcat(t8_1,t8_2),dt);
% histogram(vertcat(t8_1,t8_2),100)
title(['8 nm steps, tau=' num2str(1/fp.k1)])

fx=Hypoexponentialpdf(nt,ks(1:midp));
plot(nt,fx*N*dt*2,'r')


function [ft,nt,ht]=plotwithfit(dat,dt)
n=0:dt:max(dat);

%double exp fit:
ht=histcounts(dat,n);
nt=n(1:end-1)+dt/2;
[htmax,imax]=max(ht);
timefun=@(k1,A,x) A*k1^2*x.*exp(-k1*x);
ks=1/nt(imax);
As=1/ks*htmax*2;

startp=[ks, As];
maxtime=inf;
mintime=0;
indt=nt<maxtime;
indt=indt & nt>mintime;
ft=fit(nt(indt)',ht(indt)',timefun,'StartPoint',startp);
hold off
stairs(nt,ht,'k')
hold on
% plot(nt,timefun(startp(1),startp(2),nt))
plot(nt,ft(nt),'r--')
end
