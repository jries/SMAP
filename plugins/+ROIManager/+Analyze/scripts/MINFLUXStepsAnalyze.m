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
ds=.25;
n=round(min(stepsize)):ds:max(stepsize);
histogram(stepsize,n,'DisplayStyle','stairs')
xlabel('stepsize(nm)')
hs=histcounts(stepsize,n);
nf=n(1:end-1)+ds/2;
xlim([0 quantile(stepsize,.995)])
%Gauss fit
fs=fit(nf',hs','gauss1');
hold on
plot(nf,fs(nf),'r')
hold off
ff='%2.1f';
title(['<step>=' num2str(fs.b1,ff) ', sigma=' num2str(fs.c1,ff) ',mean(step)=' num2str(mean(stepsize),ff) ', std=' num2str(std(stepsize),ff)])

dt=1.;
n=0:dt:max(steptime);
subplot(1,2,2)
hold off
histogram(steptime,n,'DisplayStyle','stairs')
xlabel('step time (ms)')
maxtime=quantile(steptime,.98);
xlim([0 maxtime])

%double exp fit:
ht=histcounts(steptime,n);
nt=n(1:end-1)+dt/2;

timefun=@(k1,k2,A,x) A*(exp(-k1*x)-exp(-k2*x));
startp=[1/mean(steptime), 10/mean(steptime), max(ht)];
hold on
% plot(nt,timefun(startp(1),startp(2),startp(3),nt),'g')
indt=nt<maxtime;
ft=fit(nt(indt)',ht(indt)',timefun,'StartPoint',startp)
plot(nt,ft(nt),'r--')
[~,imax]=max(ht);
ftx=fit(nt(imax:end)',ht(imax:end)','exp1')
plot(nt(imax:end),ftx(nt(imax:end)),'r')
title(['1/k = ' num2str(-1/ftx.b,ff) ' ms, robustmean(steptime) = ' num2str(robustMean(steptime),ff) ' N = ' num2str(length(steptime),4)])
