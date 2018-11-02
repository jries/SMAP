% function analyzeBurstsIva2
%burst analysis for Alla
global path
%parameters to set


isczi=false; %load zeiss czi files, otherwise load txt files

exchangechannels=false; %if to exchange A and B. A needs to decay faster. false/true, 0/1
% pointspersecond=1000;%raw data points per second;
timestepsdecay=20; %how many data points to calculate to fit decays.
cutofffactor=3; %a relative factor deciding if to only take strong peaks or also weak ones.
if isczi
    selectstr='*.czi';
else
     selectstr='*.txt';
end

[file, path]=uigetfile([path filesep selectstr],'MultiSelect','on');
if ~iscell(file)
    file={file};
end

halfTimeMatrix = zeros(length(file), 3);

for fileIndex=1:length(file)
    int1=[];
    int2=[];
    
    if isczi
        rr=bfopen([path file{fileIndex}]);
        [M,meta]=readCzi(rr);

        M1=M(:,:,:,1);
        M2=M(:,:,:,2);

        int1h=M1(:);
        int2h=M2(:);
        pointspersecond=length(int1h)/meta.totaltime;
    else
        
        fid=fopen([path file{fileIndex}]);
        fgetl(fid);fgetl(fid);
        dat=textscan(fid,'%f\t%f\t%f\t%f');
        fclose(fid);
        int1h=dat{2};
        int2h=dat{4};
        pointspersecond=1/(dat{1}(2)-dat{1}(1));
    end
    
    if exchangechannels
        int1=vertcat(int1,int2h);
        int2=vertcat(int2,int1h);
    else
        int1=vertcat(int1,int1h);
        int2=vertcat(int2,int2h);
    end

    
   

    %calculate cutoff
    m1=quantile(int1(~isnan(int1)), 0.9);
    m2=quantile(int2(~isnan(int2)), 0.9);

    s1=nanstd(int1(int1<m1));
    s2=nanstd(int2(int2<m2));

    co1=nanmedian(int1)+cutofffactor*s1;
    co2=nanmedian(int2)+cutofffactor*s2;


    timewindow=round(length(int1)/timestepsdecay);
    timepoints=(1:timewindow:length(int1));
    timepointss=timepoints(1:end-1)/pointspersecond;

    %membrane: in channel 2
    ratiob=[];
    i1=[];
    i2=[];
    i1b=[];
    i2b=[];
    bg1=quantile(int1(~isnan(int1)),0.1);
    bg2=quantile(int2(~isnan(int2)),0.1);
        
    for k=length(timepoints)-1:-1:1
        timeind=false(size(int1));
        timeind(timepoints(k):timepoints(k+1))=true;
        
        int1h=int1(timeind)-bg1;
        int2h=int2(timeind)-bg2;


        inburst=int2h>co2; %only use membranes to define burst
        int1b=int1h(inburst);
        int2b=int2h(inburst);

        i1(k)=mean(int1h);
        i2(k)=mean(int2h);
        i1b(k)=mean(int1b);
        i2b(k)=mean(int2b);

        ratiob(k)=mean(int1b./int2b);

    end
    figure('Name', file{fileIndex})
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


    % half time BEGIN
    subplot(2,2,4)
    measurementtime=300;
    fitfun=@(a,b,x)(a*exp(-b*x));
    titelstr={};

    if sum(isinf(i1b)|isnan(i1b)|isinf(i2b)|isnan(i2b)) < 1
        halfTime1=fit(timepointss', (i1b./i2b)',fitfun,'StartPoint',[1, 2/measurementtime],'Lower',[0 -1e5],'Upper',[1e8 1e5 ],'Robust','LAR');
        titlestr{1}=['half time <i1b>./<i2b> = ' num2str(1/halfTime1.b, 4) 's'];
        halfTimeMatrix(fileIndex, 1) = 1/halfTime1.b;
    else
        titlestr{1}='half time <i1b>./<i2b> = error';
    end
    
    halfTime2=fit(timepointss', ratiob',    fitfun,'StartPoint',[1, 2/measurementtime],'Lower',[0 -1e5],'Upper',[1e8 1e5 ],'Robust','LAR');
    titlestr{2}=['half time <i1b/i2b>    = ' num2str(1/halfTime2.b, 4) 's'];
    halfTimeMatrix(fileIndex, 2) = 1/halfTime2.b;
    
    if sum(isinf(i1)|isnan(i1)|isinf(i2)|isnan(i2)) < 1
        halfTime3=fit(timepointss', (i1./i2)',  fitfun,'StartPoint',[1, 2/measurementtime],'Lower',[0 -1e5],'Upper',[1e8 1e5 ],'Robust','LAR');
        titlestr{3}=['half time <i1>./<i2>   = ' num2str(1/halfTime3.b, 4) 's'];
        halfTimeMatrix(fileIndex, 3) = 1/halfTime3.b;
    else
        titlestr{3}='half time <i1>./<i2>   = error';
    end

    title(titlestr);
    % half time END

    figure('Name', file{fileIndex})
    time=(1:length(int1))/pointspersecond;
    plot(time,int1-bg1);hold on
    plot(time,int2-bg2);
    plot([time(1) time(end)],[co2 co2]-bg2)


    hold off
    legend('I1','I2','Cutoff2')
end

% write to file
csvwrite("~/Desktop/test.csv", halfTimeMatrix)