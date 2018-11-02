% function analyzeBurstsIva2
%burst analysis for Iva
global path

%parameters to set
measurementtime=300; %measuremetn time in seconds to calculate the length of the data
exchangechannels=false; %if to exchange A and B. A needs to decay faster. false/true, 0/1
binning=1; %how many data points to bin for analysis
burstquantile=0.75; %fraction of data points which are not in bursts.
timestepsdecay=10; %how many data points to calculate to fit decays.
pointspersecond=1000;%raw data points per second;
timestepscatter=3; %scatter plot: how many time windows
%read the data
maxread=measurementtime*pointspersecond-100;
[file, path]=uigetfile([path filesep '*.txt']);
dat=dlmread([path file],'\t',[2,0,maxread,3]);

if exchangechannels
    int1=dat(:,4);
    int2=dat(:,2);
else
    int1=dat(:,2);
    int2=dat(:,4);
end

int1b=bintrace(int1,binning);
int2b=bintrace(int2,binning);

%find bursts
cutoff1=myquantile(int1b,burstquantile);
cutoff2=myquantile(int2b,burstquantile);




burst=int1b>cutoff1|int2b>cutoff2;

% %%
%scatter plot

dt=round(length(burst)/(timestepscatter+1));
timerange=1:dt:length(burst);
figure('Name', 'Arf intensity in burst against coatomer intensity in bursts', ...
    'NumberTitle', 'off');
hold off
plotsymbols = {'x', 'ok', 'dr', '<g'};
for k=1:length(timerange)-1
    timeind=false(size(burst));
    timeind(timerange(k):timerange(k+1))=true;
    indgood=burst&timeind;
    plot(log(int1b(indgood)),log(int2b(indgood)),plotsymbols{k});
    hold on
end
xlabel('log int2')
ylabel('log int1')
legend('t1','t2','t3','t4')
title(corr(int1b(indgood),int2b(indgood)));
%%
%time trace

figure('Name', 'Intensity trace', ...
    'NumberTitle', 'off')
plot(int1b,'-+r')
hold on
plot(int2b,'-og')
hold off

%time trace seperately
figure('Name', 'Intensity trace Arf', ...
    'NumberTitle', 'off')
plot(int1b,'-+r')
figure('Name', 'Intensity trace coatomer', ...
    'NumberTitle', 'off')
plot(int2b,'-og')


%%
%decay curves

timewindow=round(length(int1b)/timestepsdecay);
timepoints=(1:timewindow:length(int1b));
timepointss=timepoints(1:end-1)*binning/pointspersecond;

correlation=zeros(length(timepoints)-1,1);
mycorr=zeros(length(timepoints)-1,1);
rmean=zeros(length(timepoints)-1,1);
rmeanburst=zeros(length(timepoints)-1,1);
n1=zeros(length(timepoints)-1,1);
n2=zeros(length(timepoints)-1,1);
int1m=zeros(length(timepoints)-1,1);
int2m=zeros(length(timepoints)-1,1);
mycorr2=mycorr;

co1=cutoff1*4;
co2=cutoff2*4;

bg1=zeros(length(timepoints)-1,1);
bg2=zeros(length(timepoints)-1,1);
rmeanburstbgc=zeros(length(timepoints)-1,1);
quantilebackgroud=0.5; %this we can change 


for k=1:length(timepoints)-1

    timeind=false(size(burst));
    timeind(timepoints(k):timepoints(k+1))=true;
    indgood=timeind;
    int1h=int1b(indgood);
    int2h=int2b(indgood);
%     co1=myquantile(int1h,burstquantile);
%     co2=myquantile(int2h,burstquantile);
    co1a(k)=co1;
    co2a(k)=co2;
    bg1(k)=quantile(int1h,quantilebackgroud);
    bg2(k)=quantile(int2h,quantilebackgroud);
    
    co1h=bg1(k)+2.5*std(int1h<bg1(k));
    co2h=bg2(k)+2.5*std(int2h<bg2(k));
co1h=co1;co2h=co2;
    inboth=int1h>co1h&int2h>co2h;
    ina=int1h>co1h;
    inb=int2h>co2h;
    ineither=ina|inb;
    
    nboth=sum(inboth);
    na=sum(ina);
    nb=sum(inb);
     mycorr(k)=nboth/(na+nb-nboth);
     correlation(k)=corr(sqrt(int1h),sqrt(int2h));
     rhere=(int1h)./(int2h);
     int1m(k)=mean(int1h);
     int2m(k)=mean(int2h);
     int1bgc=int1h-bg1(k);
     int2bgc=int2h-bg2(k);
     rherebgc=(int1bgc)./(int2bgc);

     rmean(k)=mean(rhere(~isinf(rhere)&~isnan(rhere)));
     rmeanbgc(k)=mean(rherebgc(~isinf(rherebgc)&~isnan(rherebgc)));
     rburst=(rhere(inb));
     rburstbgc=(rherebgc(inb)); %only take those where we have intensity in the second chanenl (otherwise 1/0)
     rmeanburst(k)=mean(rburst(~isinf(rburst)&~isnan(rburst)));
%      rmeanburstbgc(k)=mean(rburstbgc(~isinf(rburstbgc)&~isnan(rburstbgc)));
%     %bgc: divsision by zero and negative valuse: ratio of averages instead
%     %of average of ratios
    rmeanburstbgc(k)=mean(int1bgc(ineither))/mean(int2bgc(ineither));
%     rmeanburst(k)=mean(int1h(ineither))/mean(int2h(ineither));
%     rmean(k)=mean(int1h)/mean(int2h);
%     
     n1(k)=na;
     n2(k)=nb;
end
figure(102);hold off
nxx=1:length(bg1);
plot(nxx,bg1);hold on; plot(nxx,bg2);plot(nxx,co1a,nxx,co2a);
plot([nxx(1) nxx(end)],[cutoff1, cutoff1],[nxx(1) nxx(end)],[cutoff2, cutoff2])
legend('bg1','bg2','co1','co2');
figure(90)
subplot(1,3,1)

minv1=myquantile(int1h,0.05);
minv2=myquantile(int2h,0.05);
minr=minv1/minv2;
hold off
plot(timepointss,rmean,'g*')
hold on
plot(timepointss,rmeanburstbgc,'r*')

% plot(timepointss,rmean-minr,'m*')
plot(timepointss,rmeanburst,'b*')

legend('intensity ratio all','intensity ratio bursts bg corr','intensity ratio bursts ')

% fitfun=@(a,b,c,x)(a*exp(-b*x)+c);
fitfunc = fittype(@ (a,b,c,x) a*exp(-b*x)+c);

fitout=fit(timepointss',rmean,fitfunc,'StartPoint',[1, .1, 2/measurementtime],'Lower',[0 0 -1e8],'Upper',[1e8 1e5 1e8],'Robust','LAR');
% plot(timepointss,fitout(timepointss),'r')

fitfun=@(a,b,x)(a*exp(-b*x));
fitout2=fit(timepointss',rmean-minr,fitfun,'StartPoint',[1, 2/measurementtime],'Lower',[0 -1e5 ],'Upper',[1e8 1e5],'Robust','LAR');
% plot(timepointss,fitout2(timepointss)+minr,'m')

fitout3=fit(timepointss',rmeanburst,fitfun,'StartPoint',[1, 2/measurementtime],'Lower',[0 -1e5],'Upper',[1e8 1e5 ],'Robust','LAR');
plot(timepointss,fitout3(timepointss),'b')

fitout4=fit(timepointss',rmeanburstbgc,fitfun,'StartPoint',[1, 2/measurementtime],'Lower',[0 -1e5],'Upper',[1e8 1e5 ],'Robust','LAR');
plot(timepointss,fitout4(timepointss),'r')

fitout5=fit(timepointss',rmeanburst,fitfunc,'StartPoint',[1, 2/measurementtime, 0],'Lower',[0 0 -inf],'Upper',[1e8 1e5 inf],'Robust','LAR');
plot(timepointss,fitout5(timepointss),'m')

xlabel('time (s)');
ylabel('intensity1 / intensity 2');
titelstr={};
% titlestr{1}=['tr all = ' num2str(1/fitout.b,4) 's'];
% titlestr{2}=['tr all = ' num2str(1/fitout2.b,4) 's'];
titlestr{1}=['tr burst = ' num2str(1/fitout3.b,4) 's'];
titlestr{2}=['tr burst bgc = ' num2str(1/fitout4.b,4) 's'];
titlestr{3}=['tr burst off = ' num2str(1/fitout5.b,4) 's'];
title(titlestr)
tout(1)=1/fitout.b;tout(2)=1/fitout2.b;tout(3)=1/fitout3.b;

subplot(1,3,2)
hold off
plot(timepointss,n1,'r*')
hold on
plot(timepointss,n2,'b*')

legend('# bursts 1','# bursts 2')
xlabel('time')
ylabel('#bursts')

fitfun=@(a,b,x)(a*exp(-b*x));
fitout1=fit(timepointss',n1,fitfun,'StartPoint',[1, 2/measurementtime],'Lower',[0 0 ],'Upper',[1e8 1e5 ],'Robust','LAR');
plot(timepointss,fitout1(timepointss),'r')

fitout2=fit(timepointss',n2,fitfun,'StartPoint',[1, 2/measurementtime],'Lower',[0 0 ],'Upper',[1e8 1e5 ],'Robust','LAR');
plot(timepointss,fitout2(timepointss),'b')

titlestr2{1}=['t1 (bursts) = ' num2str(1/fitout1.b,4)];
titlestr2{2}=['t2 (bursts) = ' num2str(1/fitout2.b,4)];
title(titlestr2)

tout(4)=1/fitout1.b;tout(5)=1/fitout2.b;

subplot(1,3,3)
hold off
plot(timepointss,int1m,'r*')
hold on

plot(timepointss,int2m,'b*')

legend('mean intensity1','mean intensity 2')
xlabel('time')
ylabel('intensity')

fitfun=@(a,b,x)(a*exp(-b*x));

fitout3=fit(timepointss',int1m,fitfun,'StartPoint',[1, 2/measurementtime],'Lower',[0 0 ],'Upper',[1e8 1e5 ],'Robust','LAR');
plot(timepointss,fitout3(timepointss),'r')
fitout4=fit(timepointss',int2m,fitfun,'StartPoint',[1, 2/measurementtime],'Lower',[0 0 ],'Upper',[1e8 1e5 ],'Robust','LAR');
plot(timepointss,fitout4(timepointss),'b')

titlestr3{1}=['t1 (intensity) = ' num2str(1/fitout3.b,4)];
titlestr3{2}=['t2 (intensity) = ' num2str(1/fitout4.b,4)];
tout(6)=1/fitout3.b;tout(7)=1/fitout4.b;
% titlestr{3}=['kr burst = ' num2str(1/fitout.b)];
title(titlestr3)

allts=sprintf('%4.2f \t',tout)
