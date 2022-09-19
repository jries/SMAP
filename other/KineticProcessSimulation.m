%Simulate kinetic processes, chain of events with kinetic constants ki

N=100000; %number of data points
ks=1./[30 10 30 10]; % list of kinetic constants
dt=0.005;

timet=zeros(N,length(ks));

for k=1:length(ks)
    timet(:,k)=exprnd(ks(k),N,1);
end

%total cycle time:
t2=sum(timet,2);
figure(88);
subplot(2,1,1)
fp=plotwithfit(t2,dt);
title(['16 nm steps, k=' num2str(fp.k1)])

subplot(2,1,2)
l=size(timet,2); 
midp=round(l/2);
t8_1=sum(timet(:,1:midp),2);
t8_2=sum(timet(:,end-midp+1:end),2);
fp=plotwithfit(vertcat(t8_1,t8_2),dt);
% histogram(vertcat(t8_1,t8_2),100)
title(['8 nm steps, k=' num2str(fp.k1)])

function ft=plotwithfit(dat,dt)
n=0:dt:max(dat);

%double exp fit:
ht=histcounts(dat,n);
nt=n(1:end-1)+dt/2;
[~,imax]=max(ht);
timefun=@(k1,A,x) A*k1^2*x.*exp(-k1*x);
startp=[1/mean(dat), 3*max(ht)*mean(dat)^2/ht(imax)];
maxtime=inf;
mintime=0;
indt=nt<maxtime;
indt=indt & nt>mintime;
ft=fit(nt(indt)',ht(indt)',timefun,'StartPoint',startp);
hold off
stairs(nt,ht)
hold on
plot(nt,ft(nt),'r--')
end
