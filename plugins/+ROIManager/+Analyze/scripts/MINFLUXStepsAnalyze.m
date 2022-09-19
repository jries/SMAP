%analysis of MINFLUX tracks
sites=g.locData.SE.sites;

stepsize=[]; %define the variables you want to read out
steptime=[];

pluginname='StepsMINFLUX';
fieldname1='steps';

vel=[];
tracklengthlocs=[];

for k=1:length(sites)
    if ~isfield(sites(k).evaluation,pluginname) || ~isfield(sites(k).evaluation.(pluginname),fieldname1) % if no steps are found, look at next site
        k
        continue
    end
    sh=sites(k).evaluation.(pluginname).(fieldname1);
    stepsize(end+1:end+length(sh.stepsize))=sh.stepsize; %add values from current site to the list
    steptime(end+1:end+length(sh.dwelltime))=sh.dwelltime;
    vel(end+1)=sites(k).evaluation.(pluginname).stattrack.velocity;
    tracklengthlocs(end+1)=sites(k).evaluation.(pluginname).statall.nlocs;
end

figure(188) %do the plotting
subplot(2,2,1)
ds=.5;
n=round(min(stepsize)):ds:max(stepsize);
hs=histcounts(stepsize,n);
stairs(n,[hs,0])
% histogram(stepsize,n,'DisplayStyle','stairs')
xlabel('stepsize(nm)')

nf=n(1:end-1)+ds/2;
mn=max(4*round(mean(stepsize)/2), max(stepsize));

xlim([0 mn])
%Gauss fit
fs=fit(nf',hs','gauss1');
hold on
np=0:ds/5:mn;
plot(np,fs(np),'r')
hold off
ff='%2.1f';
title(['fit: <step>=' num2str(fs.b1,ff) ', sigma=' num2str(fs.c1,ff) ',mean(step)=' num2str(mean(stepsize),ff) ', std=' num2str(std(stepsize),ff)])

dt=2;

n=0:dt:max(steptime);

%double exp fit:
ht=histcounts(steptime,n);
nt=n(1:end-1)+dt/2;


subplot(2,2,2)
hold off
stairs(n,[ht,0],'k')
xlabel('step time (ms)')
maxtime=quantile(steptime,.98);
xlim([0 maxtime])
hold on
%startp
[htmax,imax]=max(ht);
kstart=1/nt(imax);
indt=nt<maxtime;
indt=indt & nt>0;
Ns=sum(ht(indt));



%Hypoexp fit 16 nm
fitf=@(k1,k2,A,x) A*Ns*dt*Hypoexp4fit(x,k1,k2);
fthyp=fit(nt(indt)',ht(indt)',fitf,'StartPoint',[kstart,kstart*4,1]);
plot(nt,fthyp(nt),'b')

%Hypoexp fit 8 nm
fitf=@(k1,k2,A,x) A*Ns*dt*Hypoexponentialpdf(x,[k1,k2]);
[htmax,imax]=max(ht);
kstart=1/nt(imax);
% plot(nt,fitf(kstart,kstart*8,1,nt),'m--')
fthyp2=fit(nt(indt)',ht(indt)',fitf,'StartPoint',[kstart,kstart*8,1]);
plot(nt,fthyp2(nt),'m')

%convolution of two exponentials with same rates, Peng dynein paper
timefun=@(k1,A,x) A*Ns*dt*k1^2*x.*exp(-k1*x);
startp=[kstart,1];

ft=fit(nt(indt)',ht(indt)',timefun,'StartPoint',startp);

plot(nt,ft(nt),'g--')

%exp
ftx=fit(nt(imax:find(indt,1,'last'))',ht(imax:find(indt,1,'last'))','exp1');
plot(nt(imax:find(indt,1,'last')),ftx(nt(imax:find(indt,1,'last'))),'r')

title(['Peng: 1/k = ' num2str(1/ft.k1,ff) ' ms, exp: 1/k = ' num2str(-1/ftx.b,ff)...
    '; 4exp: ' num2str(1/fthyp.k1,ff) ',' num2str(1/fthyp.k2,ff) ...
     '; 2exp: ' num2str(1/fthyp2.k1,ff) ',' num2str(1/fthyp2.k2,ff) ...
     '; N = ' num2str(length(steptime),4)]);

legend('data',...
    ['4exp: ' num2str(1/fthyp.k1,ff) ',' num2str(1/fthyp.k2,ff)],...
    ['2exp: ' num2str(1/fthyp2.k1,ff) ',' num2str(1/fthyp2.k2,ff)],...
    ['Peng: 1/k = ' num2str(1/ft.k1,ff)],...
    ['exp: 1/k = ' num2str(-1/ftx.b,ff)])



figure(180)
subplot(1,2,1)
n=0:0.05:max(vel);
histogram(vel,n)
xlabel('velocity (nm/s)')
ff2='%2.0f';
title(['v = ' num2str(mean(vel)*1000,ff2) '±' num2str(std(vel)*1000,ff2) ' nm/s '])

subplot(1,2,2)
hold off
nl=0:500:max(tracklengthlocs);
histogram(tracklengthlocs,nl)
h=histcounts(tracklengthlocs,nl);
fp=fit(nl(1:end-1)',h','exp1');
hold on
plot(nl,fp(nl))

title(['<nlocs> = ' num2str(mean(tracklengthlocs),ff2) ', median = ' num2str(median(tracklengthlocs),ff2) ', exp = ' num2str(-1/fp.b,ff2)])


figure(188)
subplot(2,2,3)
hold off
htc=cumsum(ht);
stairs(nt,htc,'k')
xlabel('step time (ms)')
maxtime=quantile(steptime,.98);
% maxtime=100;
xlim([0 maxtime])
hold on
%startp
[htmax,imax]=max(ht);
kstart=1/nt(imax);
indt=nt<maxtime;
indt=indt & nt>0;
Ns=sum(ht(indt));



%Hypoexp fit 16 nm
fitf=@(k1,k2,A,x) A*Ns*dt*cumsum(Hypoexp4fit(x,k1,k2));
fthyp=fit(nt(indt)',htc(indt)',fitf,'StartPoint',[kstart,kstart*4,1]);
plot(nt,fthyp(nt),'b')

%Hypoexp fit 8 nm
fitf=@(k1,k2,A,x) A*Ns*dt*cumsum(Hypoexponentialpdf(x,[k1,k2]));
[htmax,imax]=max(ht);
kstart=1/nt(imax);
% plot(nt,fitf(kstart,kstart*8,1,nt),'m--')
fthyp2=fit(nt(indt)',htc(indt)',fitf,'StartPoint',[kstart,kstart*8,1]);
plot(nt,fthyp2(nt),'m')


%convolution of two exponentials with same rates, Peng dynein paper
timefun=@(k1,A,x) A*Ns*dt*k1^2*cumsum(x.*exp(-k1*x));
startp=[kstart,1];

ft=fit(nt(indt)',htc(indt)',timefun,'StartPoint',startp);

plot(nt,ft(nt),'g--')

%exp
ftx=fit(nt(imax:find(indt,1,'last'))',ht(imax:find(indt,1,'last'))','exp1');

plot(nt(imax:find(indt,1,'last')),sum(ht(1:imax-1))+cumsum(ftx(nt(imax:find(indt,1,'last')))),'r')

% title(['Peng: 1/k = ' num2str(1/ft.k1,ff) ' ms, exp: 1/k = ' num2str(-1/ftx.b,ff)...
%     '; 4exp: ' num2str(1/fthyp.k1,ff) ',' num2str(1/fthyp.k2,ff) ...
%      '; 2exp: ' num2str(1/fthyp2.k1,ff) ',' num2str(1/fthyp2.k2,ff) ...
%      '; N = ' num2str(length(steptime),4)]);

legend('data',...
    ['4exp: ' num2str(1/fthyp.k1,ff) ',' num2str(1/fthyp.k2,ff)],...
    ['2exp: ' num2str(1/fthyp2.k1,ff) ',' num2str(1/fthyp2.k2,ff)],...
    ['Peng: 1/k = ' num2str(1/ft.k1,ff)],...
    ['exp: 1/k = ' num2str(-1/ftx.b,ff)],'Location','southeast')

title('Cumulative Probability Distribution')

figure(188)
subplot(2,2,4)
hold off

stairs(nt,htc*0,'k')
xlabel('step time (ms)')
xlim([0 maxtime])
hold on

%Hypoexp fit 16 nm
plot(nt,htc'-fthyp(nt),'b')


%Hypoexp fit 8 nm
plot(nt,htc'-fthyp2(nt),'m')

%convolution of two exponentials with same rates, Peng dynein paper
plot(nt,htc'-ft(nt),'g--')

plot(nt(imax:find(indt,1,'last')),htc(imax:find(indt,1,'last'))'-sum(ht(1:imax-1))-cumsum(ftx(nt(imax:find(indt,1,'last')))),'r')

legend('data',...
    ['4exp: ' num2str(1/fthyp.k1,ff) ',' num2str(1/fthyp.k2,ff)],...
    ['2exp: ' num2str(1/fthyp2.k1,ff) ',' num2str(1/fthyp2.k2,ff)],...
    ['Peng: 1/k = ' num2str(1/ft.k1,ff)],...
    ['exp: 1/k = ' num2str(-1/ftx.b,ff)],'Location','southeast')

title('residuals CPD')
