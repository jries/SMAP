%analysis of MINFLUX tracks
sites=g.locData.SE.sites;

stepsize=[]; %define the variables you want to read out
steptime=[];

pluginname='StepsMINFLUX';
fieldname1='steps';


for k=1:length(sites)
    if ~isfield(sites(k).evaluation,pluginname) || ~isfield(sites(k).evaluation.(pluginname),fieldname1) % if no steps are found, look at next site
        continue
    end
    sh=sites(k).evaluation.(pluginname).(fieldname1);
    stepsize(end+1:end+length(sh.stepsize))=sh.stepsize; %add values from current site to the list
    steptime(end+1:end+length(sh.dwelltime))=sh.dwelltime;
end

figure(88) %do the plotting
subplot(1,2,1)
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

dt=10;
n=0:dt:max(steptime);

%double exp fit:
ht=histcounts(steptime,n);
nt=n(1:end-1)+dt/2;


subplot(1,2,2)
hold off
stairs(n,[ht,0])
% histogram(steptime,n,'DisplayStyle','stairs')
xlabel('step time (ms)')
maxtime=quantile(steptime,.98);
xlim([0 maxtime])


%two exponentials
timefun=@(k1,k2,A,x) A*(exp(-k1*x)-exp(-k2*x));
startp=[1/mean(steptime), 10/mean(steptime), max(ht)];
indt=nt<maxtime;
ft=fit(nt(indt)',ht(indt)',timefun,'StartPoint',startp);
[~,imax]=max(ht);
%convolution of two exponentials with same rates, Peng dynein paper
timefun=@(k1,A,x) A*k1^2*x.*exp(-k1*x);
startp=[2/mean(steptime), max(ht)*mean(steptime)^2/ht(imax)];
indt=nt<maxtime;
ft=fit(nt(indt)',ht(indt)',timefun,'StartPoint',startp);
hold on
plot(nt,ft(nt),'r--')

%convolution of two exponentials with different rates, 
timefun=@(k1,k2,A,x) A*(k1*k2/((k2-k1)))^2*((x-(2/((k2-k1)))).*exp(-k1*x)+(x-(2/((k2-k1)))).*exp(-k2*x));
startp=[ft.k1,ft.k1/10, max(ht)*mean(steptime)^2/ht(imax)];
ft2=fit(nt(indt)',ht(indt)',timefun,'StartPoint',startp);
hold on
plot(nt,ft2(nt),'g--')

ftx=fit(nt(imax:find(indt,1,'last'))',ht(imax:find(indt,1,'last'))','exp1');
plot(nt(imax:find(indt,1,'last')),ftx(nt(imax:find(indt,1,'last'))),'r')
title(['Peng: 1/k = ' num2str(1/ft.k1,ff) ' ms, exp: 1/k = ' num2str(-1/ftx.b,ff) ' N = ' num2str(length(steptime),4) ...
    ', 2exp: 1/k1,: ' num2str(1./ft2.k1,ff) ',1/k2,: ' num2str(1./ft2.k2,ff)])
% title(['Peng: 1/k = ' num2str(-1/ft.k1,ff) ' ms, robustmean(steptime) = ' num2str(robustMean(steptime),ff) ' N = ' num2str(length(steptime),4)])


