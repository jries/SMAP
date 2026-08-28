%analysis of MINFLUX tracks

% Adjust the histogram time bin, dt, for best dwell time fit.
dt=4;

% groupfield='groupindex' or 'tid'; %Choose either of the two depending on how you grouped locs into tracks
groupfield='tid';

sites=g.locData.SE.sites;

stepsize=[]; %define the variables you want to read out
steptime=[];

pluginname='StepsMINFLUX';
fieldname1='steps';
fieldname2='statall';

vel=[];
tracklengthlocs=[];
badsteps=[];
indstep=[];
newtrack=[];
dtmedian = [];
velocity = [];
nlocs = [];
tracklength = [];
stddetrendx = [];
stddetrendy = [];
stddetrendz = [];
msdall= [];
dwelltime = {};
timeall = {};
deltaT = {};
efoall = {};
cfrall = {};


for k=1:length(sites)
    if ~isfield(sites(k).evaluation,pluginname) || ~isfield(sites(k).evaluation.(pluginname),fieldname1) % if no steps are found, look at next site
        k
        continue
    end
    if sites(1,k).annotation.use == 1
    sh=sites(k).evaluation.(pluginname).(fieldname1);
    sh2=sites(k).evaluation.(pluginname).(fieldname2);
    
    stepsize(end+1:end+length(sh.stepsize))=sh.stepsize; %add values from current site to the list
    steptime(end+1:end+length(sh.dwelltime))=sh.dwelltime;
    stddetrendx(end+1:end+length(sh.stddetrend.x))=sh.stddetrend.x;
    stddetrendy(end+1:end+length(sh.stddetrend.y))=sh.stddetrend.y;
    if isfield(sh.stddetrend, 'z')
        stddetrendz(end+1:end+length(sh.stddetrend.z))=sh.stddetrend.z; 
    end   
    dtmedian(end+1:end+length(sh2.dtmedian))=sh2.dtmedian;
    velocity(end+1:end+length(sh2.velocity))=sh2.velocity;
    nlocs(end+1:end+length(sh2.nlocs))=sh2.nlocs;
    tracklength(end+1:end+length(sh2.tracklength))=sh2.tracklength;
    vel(end+1)=sites(k).evaluation.(pluginname).stattrack.velocity;
    tracklengthlocs(end+1)=sites(k).evaluation.(pluginname).statall.nlocs;
    id=g.locData.SE.sites(k).evaluation.StepsMINFLUX.statall.id;

    index=g.locData.loc.(groupfield)==id;

    dwelltime{k} = sh.dwelltime;
    timeall{k} = g.locData.loc.time(index);
    deltaT{k} = diff(timeall{k});
    efoall{k} = g.locData.loc.efo(index);
    cfrall{k} = g.locData.loc.cfr(index);
    
    if isfield(sh,'badsteps')
    badsteps(end+1:end+length(sh.badsteps))=sh.badsteps;
    indstep(end+1:end+length(sh.indstepglobal))=sh.indstepglobal;
    nt=0*sh.indstepglobal;
    nt(1)=1;
    nt(end)=1;
    newtrack(end+1:end+length(sh.indstepglobal))=nt;
    end
    end
   
end


loctime = [];
for kk = 1:length(deltaT)
    loctime = cat(1,loctime,double(deltaT{kk}));
end
cfrarray = [];
for kk = 1:length(cfrall)
    cfrarray = cat(1,cfrarray,double(cfrall{kk}));
end
efoarray = [];
for kk = 1:length(efoall)
    efoarray = cat(1,efoarray,double(efoall{kk}));
end

% Saving statistical data to a clipboard of the system. Simply paste to e.g. Excel sheet.
clipboard('copy', [sprintf([num2str(mean(stddetrendx, 'omitnan')) '\t' num2str(std(stddetrendx, 'omitnan')) ...
    '\t' num2str(std(stddetrendx, 'omitnan')/sqrt(length(stddetrendx))) '\t' num2str(mean(stddetrendy, 'omitnan')) ...
    '\t' num2str(std(stddetrendy, 'omitnan')) '\t' num2str(std(stddetrendy, 'omitnan')/sqrt(length(stddetrendy))) ...
    '\t' num2str(mean(stddetrendz, 'omitnan')) '\t' num2str(std(stddetrendz, 'omitnan'))...
    '\t' num2str(std(stddetrendz, 'omitnan')/sqrt(length(stddetrendz))) '\t' num2str(mean(tracklength)) '\t' num2str(std(tracklength)) ...
    '\t' num2str(mean(vel)) '\t' num2str(std(vel))  '\t' num2str(median(loctime))  '\t' num2str(median(cfrarray(cfrarray>0))) '\t' num2str(median(efoarray(efoarray>0)))  ])])

display(sprintf('Precision x mean \t std \t sem \t Precision y mean \t std \t sem \t Precision z mean \t std \t sem \t Track length mean \t std \t Average walking speed \t std \t Median localization time (ms) \t Median CFR \t Median EFO'))


%identify bad steps
% 1. bad time resolution
% 2. optionally: remove first step per track (it goes to boundary,
% not necessarily a real step, last step already not taken into account)
%instep refers to end of step
removefirst=0;
removemintime=1;

window=2; %+/- window
mintime=4; %ms

badind=false(size(indstep));
for k=1:length(indstep)
    if newtrack(k) && removefirst
        badind(k)=true;
    else
        tstep=g.locData.loc.time(indstep(k)-window:indstep(k)+window);
        if removemintime && any(diff(tstep)>mintime)
            badind(k)=true;
            badind(k+1)=true;
        end
    end
end

% stepsize(badind)=[];
steptime(badind(1:length(indstep)))=[];

figure(188) %do the plotting
f=gcf;f.Renderer='painters';
subplot(2,3,1)
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
ffv='%2.0f';
title(['fit: <step>=' num2str(fs.b1,ff) ', sigma=' num2str(fs.c1,ff) ',mean(step)=' num2str(mean(stepsize),ff) ', std=' num2str(std(stepsize),ff)])
stepsizemean=fs.b1;


n=0:dt:max(steptime);

%double exp fit:
ht=histcounts(steptime,n);
nt=n(1:end-1)+dt/2;


subplot(2,3,2)
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
% fitf=@(k1,k2,A,x) A*Ns*dt*Hypoexp4fit(x,k1,k2);
fitf=@(k1,k2,A,x) A*Ns*dt*Hypoexponentialpdf(x,[k1,k2,k1,k2]);
fthyp=fit(nt(indt)',ht(indt)',fitf,'StartPoint',[kstart,kstart*18,1]);
plot(nt,fthyp(nt),'r')
meandtfit4e2=sum(fthyp(nt).*nt')/sum(fthyp(nt));
v4e2=stepsizemean/meandtfit4e2*1000; %um/ms

%16 nm, four equal constants
% fitf=@(k1,A,x) A*Ns*dt*Hypoexp4fit(x,k1,k1);
fitf=@(k1,A,x) A*Ns*dt*Erlangdistribution(x,k1,4);
fthyp4_1=fit(nt(indt)',ht(indt)',fitf,'StartPoint',[kstart,1]);
plot(nt,fthyp4_1(nt),'b-.')
meandtfit4e1=sum(fthyp4_1(nt).*nt')/sum(fthyp4_1(nt));
v4e1=stepsizemean/meandtfit4e1*1000; %um/ms

% 
% %convolution of two exponentials with same rates, Peng dynein paper
% % timefun=@(k1,A,x) A*Ns*dt*k1^2*x.*exp(-k1*x);
% timefun=@(k1,A,x) A*Ns*dt*Erlangdistribution(x,k1,2);
% startp=[kstart,1];
% ft=fit(nt(indt)',ht(indt)',timefun,'StartPoint',startp);
% plot(nt,ft(nt),'b--')
% meandtfit2e1=sum(ft(nt).*nt')/sum(ft(nt));
% v2e1=stepsizemean/meandtfit2e1*1000; %um/ms

%Hypoexp fit 8 nm
fitf=@(k1,k2,A,x) A*Ns*dt*Hypoexponentialpdf(x,[k1,k2]);
[htmax,imax]=max(ht);
kstart=1/nt(imax);
% plot(nt,fitf(kstart,kstart*8,1,nt),'m--')
fthyp2=fit(nt(indt)',ht(indt)',fitf,'StartPoint',[kstart/2,kstart*2,1]);
% fthyp2=fit(nt(indt)',ht(indt)',fitf,'StartPoint',[ft.k1*.9,ft.k1*1.1,1]);
plot(nt,fthyp2(nt),'m')
meandtfit2e2=sum(fthyp2(nt).*nt')/sum(fthyp2(nt));
v2e2=stepsizemean/meandtfit2e2*1000; %um/ms

% %3 exp
% timefun=@(k1,A,x) A*Ns*dt*Erlangdistribution(x,k1,3);
% startp=[kstart,1];
% ft3=fit(nt(indt)',ht(indt)',timefun,'StartPoint',startp);
% plot(nt,ft3(nt),'b:')
% meandtfit3e1=sum(ft3(nt).*nt')/sum(ft3(nt));
% v3e1=stepsizemean/meandtfit3e1*1000; %um/ms

%exp
ftx=fit(nt(imax:find(indt,1,'last'))',ht(imax:find(indt,1,'last'))','exp1');
plot(nt(imax:find(indt,1,'last')),ftx(nt(imax:find(indt,1,'last'))),'c')
meandtfit1e1=sum(ftx(nt).*nt')/sum(ftx(nt));
v1e1=stepsizemean/meandtfit1e1*1000; %um/ms

title(['N = ' num2str(length(steptime),4)]);

legend('data, shown is 1/k, v: nm/s, dt: ms',...
    ['4exp 2 k: ' num2str(1/fthyp.k1,ff) ',' num2str(1/fthyp.k2,ff) ', v=' num2str(v4e2,ffv) ', dt=' num2str(meandtfit4e2,ff)],...
    ['4exp 1 k: ' num2str(1/fthyp4_1.k1,ff) ', v=' num2str(v4e1,ffv) ', dt=' num2str(meandtfit4e1,ff)],...
    ['2exp: ' num2str(1/fthyp2.k1,ff) ',' num2str(1/fthyp2.k2,ff) ', v=' num2str(v2e2,ffv) ', dt=' num2str(meandtfit2e2,ff)],...
    ['1exp: ' num2str(-1/ftx.b,ff) ', v=' num2str(v1e1,ffv) ', dt=' num2str(meandtfit1e1,ff)])

%     ['2exp 1 k: ' num2str(1/ft.k1,ff) ', v=' num2str(v2e1,ffv) ', dt=' num2str(meandtfit2e1,ff)],...
%     ['3exp 1 k: ' num2str(1/ft3.k1,ff) ', v=' num2str(v3e1,ffv) ', dt=' num2str(meandtfit3e1,ff)],...

figure(188)
subplot(2,6,11)
% nv=0:0.05:max(vel);
nv=0:max(vel)/20:max(vel);
histogram(vel,nv)
xlabel('velocity (nm/s)')
ff2='%2.0f';
title(['v = ' num2str(mean(vel)*1000,ff2) '±' num2str(std(vel)*1000,ff2) ' nm/s '])

subplot(2,6,12)
hold off
% nl=0:500:max(tracklengthlocs);
nl=0:max(tracklengthlocs)/20:max(tracklengthlocs);
histogram(tracklengthlocs,nl)
h=histcounts(tracklengthlocs,nl);
fp=fit(nl(1:end-1)',h','exp1');
hold on
plot(nl,fp(nl))

title(['<nlocs> = ' num2str(mean(tracklengthlocs),ff2) ', median = ' num2str(median(tracklengthlocs),ff2) ', exp = ' num2str(-1/fp.b,ff2)])


figure(188)
subplot(2,3,3)
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
fitf=@(k1,k2,A,x) A*Ns*dt*cumsum(Hypoexponentialpdf(x,[k1,k2,k1,k2]));
fthypc=fit(nt(indt)',htc(indt)',fitf,'StartPoint',[kstart,kstart*4,1]);
plot(nt,fthypc(nt),'b')

% fitf=@(k1,A,x) A*Ns*dt*cumsum(Hypoexp4fit(x,k1,k1));
fitf=@(k1,A,x) A*Ns*dt*cumsum(Erlangdistribution(x,k1,4));
fthypc4_1=fit(nt(indt)',htc(indt)',fitf,'StartPoint',[kstart,1]);
plot(nt,fthypc4_1(nt),'c--')

%convolution of two exponentials with same rates, Peng dynein paper
% timefun=@(k1,A,x) A*Ns*dt*k1^2*cumsum(x.*exp(-k1*x));
timefun=@(k1,A,x) A*Ns*dt*cumsum(Erlangdistribution(x,k1,2));
startp=[kstart,1];
ft=fit(nt(indt)',htc(indt)',timefun,'StartPoint',startp);
plot(nt,ft(nt),'g')


%Hypoexp fit 8 nm
fitf=@(k1,k2,A,x) A*Ns*dt*cumsum(Hypoexponentialpdf(x,[k1,k2]));
[htmax,imax]=max(ht);
kstart=1/nt(imax);
% plot(nt,fitf(kstart,kstart*8,1,nt),'m--')
fthyp2c=fit(nt(indt)',htc(indt)',fitf,'StartPoint',[ft.k1*0.9,ft.k1*1.1,1]);
plot(nt,fthyp2c(nt),'m--')




%exp
ftx=fit(nt(imax:find(indt,1,'last'))',ht(imax:find(indt,1,'last'))','exp1');

plot(nt(imax:find(indt,1,'last')),sum(ht(1:imax-1))+cumsum(ftx(nt(imax:find(indt,1,'last')))),'r')

% title(['Peng: 1/k = ' num2str(1/ft.k1,ff) ' ms, exp: 1/k = ' num2str(-1/ftx.b,ff)...
%     '; 4exp: ' num2str(1/fthyp.k1,ff) ',' num2str(1/fthyp.k2,ff) ...
%      '; 2exp: ' num2str(1/fthyp2.k1,ff) ',' num2str(1/fthyp2.k2,ff) ...
%      '; N = ' num2str(length(steptime),4)]);

legend('data',...
    ['4exp 2 k: ' num2str(1/fthyp.k1,ff) ',' num2str(1/fthypc.k2,ff)],...
    ['4exp 1 k: ' num2str(1/fthypc4_1.k1,ff)],...
    ['Peng: 1/k = ' num2str(1/ft.k1,ff)],...
    ['2exp: ' num2str(1/fthyp2.k1,ff) ',' num2str(1/fthyp2.k2,ff)],...
    ['exp: 1/k = ' num2str(-1/ftx.b,ff)],'Location','southeast')

title('Cumulative Probability Distribution')

figure(188)
subplot(2,3,4)
hold off

stairs(nt,htc*0,'k')
xlabel('step time (ms)')
xlim([0 maxtime])
hold on

%Hypoexp fit 16 nm
plot(nt,htc'-fthypc(nt),'b')

plot(nt,htc'-fthypc4_1(nt),'c--')

%convolution of two exponentials with same rates, Peng dynein paper
plot(nt,htc'-ft(nt),'g')

%Hypoexp fit 8 nm
plot(nt,htc'-fthyp2c(nt),'m--')



plot(nt(imax:find(indt,1,'last')),htc(imax:find(indt,1,'last'))'-sum(ht(1:imax-1))-cumsum(ftx(nt(imax:find(indt,1,'last')))),'r')

legend('data',...
    ['4exp 2 k: ' num2str(1/fthyp.k1,ff) ',' num2str(1/fthypc.k2,ff)],...
    ['4exp 1 k: ' num2str(1/fthypc4_1.k1,ff)],...
     ['Peng: 1/k = ' num2str(1/ft.k1,ff)],...
    ['2exp: ' num2str(1/fthyp2.k1,ff) ',' num2str(1/fthyp2.k2,ff)],...
    ['exp: 1/k = ' num2str(-1/ftx.b,ff)],'Location','southeast')

title('residuals CPD')


subplot(2,3,5)

hold off
stairs(n,[ht,0],'k')
xlabel('step time (ms)')
maxtime=quantile(steptime,.98);
xlim([0 maxtime])
hold on


%Hypoexp fit 16 nm
plot(nt,fthypc.A*Ns*dt*Hypoexponentialpdf(nt,[fthypc.k1,fthypc.k2,fthypc.k1,fthypc.k2]),'r')


%16 nm, four equal constants
plot(nt,fthypc4_1.A*Ns*dt*Erlangdistribution(nt,fthypc4_1.k1,4),'b-.')


% %convolution of two exponentials with same rates, Peng dynein paper
% % timefun=@(k1,A,x) A*Ns*dt*k1^2*x.*exp(-k1*x);
% timefun=@(k1,A,x) A*Ns*dt*Erlangdistribution(x,k1,2);
% startp=[kstart,1];
% ft=fit(nt(indt)',ht(indt)',timefun,'StartPoint',startp);
% plot(nt,ft(nt),'b--')
% meandtfit2e1=sum(ft(nt).*nt')/sum(ft(nt));
% v2e1=stepsizemean/meandtfit2e1*1000; %um/ms

%Hypoexp fit 8 nm
plot(nt,fthyp2.A*Ns*dt*Hypoexponentialpdf(nt,[fthyp2.k1,fthyp2.k2]),'m')
title('histograms from cumsum')

legend('data',...
    ['4exp 2 k: ' num2str(1/fthyp.k1,ff) ',' num2str(1/fthypc.k2,ff)],...
    ['4exp 1 k: ' num2str(1/fthypc4_1.k1,ff)],...
    ['2exp: ' num2str(1/fthyp2.k1,ff) ',' num2str(1/fthyp2.k2,ff)])

