% name='room1';
name=inputdlg('name / room');
name=name{1};
path='/Volumes/t2ries/projects/4Pi/201905_vibrations/';

%record and show power spectrum
updatetime=10; %s =1/dF

numberOfBlocks=10;

fmax=105; %Hz for plotting only
fmin=30;
nbits=16;%
Fs=5000; %sampling frequency

%s
rec=audiorecorder(Fs,nbits,1);

figure(88);
ax=gca;
blockav=zeros(Fs*updatetime/2+1,1);
alltr=zeros(Fs*updatetime/2+1,numberOfBlocks);
for k=1:numberOfBlocks
    rec.recordblocking(updatetime);
    a=rec.getaudiodata;
    L=length(a);
    spectrum=abs(fft(a)/L);
    f=Fs*(0:(L/2))/L;
    sp1=spectrum(1:L/2+1);
    blockav=blockav+sp1;
    alltr(:,k)=sp1;
    plot(ax,f,sp1,f,blockav/k)
    legend(ax,'block','average')
    xlim([fmin fmax]);
end
averagef=blockav/k;
frequency=f;
allf=alltr;

save([path name '.mat'], 'frequency','allf','averagef');
disp('done')