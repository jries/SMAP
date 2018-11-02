% function analyzeBurstsIva2
%burst analysis for Alla
global path
%parameters to set
exchangechannels=false; %if to exchange A and B. A needs to decay faster. false/true, 0/1
% pointspersecond=1000;%raw data points per second;
timestepsdecay=20; %how many data points to calculate to fit decays.
cutofffactor=3; %a relative factor deciding if to only take strong peaks or also weak ones.

[file, path]=uigetfile([path filesep '*.txt'],'MultiSelect','on');
if ~iscell(file)
    file={file};
end

int1=[];
int2=[];
for k=1:length(file)
    fid=fopen([path file{k}]);
    fgetl(fid);fgetl(fid);
    dat=textscan(fid,'%f\t%f\t%f\t%f');
    if exchangechannels
        int1=vertcat(int1,dat{4});
        int2=vertcat(int2,dat{2});
    else
        int1=vertcat(int1,dat{2});
        int2=vertcat(int2,dat{4});
    end
end
fclose(fid);

pointspersecond=1/(dat{1}(2)-dat{1}(1));

%calculate cutoff
m1=quantile(int1(~isnan(int1)),.9);
m2=quantile(int2(~isnan(int2)),.9);
s1=nanstd(int1(int1<m1));
s2=nanstd(int2(int2<m2));

co1=nanmedian(int1)+cutofffactor*s1;
co2=nanmedian(int2)+cutofffactor*s2;


timewindow=round(length(int1)/timestepsdecay);
timepoints=(1:timewindow:length(int1));
timepointss=timepoints(1:end-1)/pointspersecond;

%membrane: in channel 2
ratiob=[];i1=[];i2=[];i1b=[];i2b=[];
for k=length(timepoints)-1:-1:1
    timeind=false(size(int1));
    timeind(timepoints(k):timepoints(k+1))=true;
    int1h=int1(timeind);int2h=int2(timeind);
    inburst=int2h>co2; %only use membranes to define burst
    int1b=int1h(inburst);int2b=int2h(inburst);
    i1(k)=mean(int1h);i2(k)=mean(int2h);
    i1b(k)=mean(int1b);i2b(k)=mean(int2b);
    ratiob(k)=mean(int1b./int2b);
    
end
figure(33);
subplot(2,2,1)
plot(timepointss,i1,timepointss,i2)
title('intensity')
legend('I1','I2')
subplot(2,2,2)
plot(timepointss,i1b,timepointss,i2b)
title('intensity in burst')
legend('I1','I2')
subplot(2,2,3)
plot(timepointss, i1b./i2b,timepointss,ratiob,timepointss,i1./i2)
title('intensity ratio')
legend('!<I1b>/<I2b>','<I1b/I2b>','<I1>/<I2>')
figure(34)
time=(1:length(int1))/pointspersecond;
plot(time,int1);hold on
plot(time,int2);
plot([time(1) time(end)],[co2 co2])
hold off
legend('I1','I2','Cutoff2')
