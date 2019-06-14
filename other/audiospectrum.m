% set audiorecorder to external device if existent
    audinf = audiodevinfo;
    if size(audinf.input, 2) > 1
        audID = audinf.input(2).ID; % get ID from external audio device
        disp(['Recording from audio device ' audinf.input(2).Name '.'])
    else
        audID = audinf.input(1).ID; % get ID from internal audio device
        disp(['Recording from audio device ' audinf.input(1).Name '.'])
    end % if
    

% name='room1';
name=inputdlg('name / room');
name=name{1};
path='/Users/diekmann/Documents/24 Audio spectra/';

%record and show power spectrum
updatetime=10; %s =1/dF

numberOfBlocks=10;

fmax=105; %Hz for plotting only
fmin=30;
nbits=16;%
Fs=5000; %sampling frequency

try
    rec=audiorecorder(Fs,nbits,1,audID);

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

        disp(['Progress: ' num2str(k*updatetime) 's of ' num2str(numberOfBlocks*updatetime) 's'])
    end

    figure
    subplot(2,1,1)
        plot(f,blockav/k, 'Color', [0.8500, 0.3250, 0.0980])
        legend('average')
        xlim([fmin fmax]);
        xlabel('Frequency (Hz)')
        ylabel('Amplitude (a.u.)')
        title(name)
    subplot(2,1,2)
        loglog(f,blockav/k, 'Color', [0.8500, 0.3250, 0.0980])
        legend('average')
        %xlim([fmin fmax]);
        xlabel('Frequency (Hz)')
        ylabel('Amplitude (a.u.)')
        title(name)   
    savefig([path 'spec__' name '.fig'])


    averagef=blockav/k;
    frequency=f;
    allf=alltr;

    save([path 'spec__' name '.mat'], 'frequency','allf','averagef');
    disp('Done. Results saved')
catch
    averagef=blockav/k;
    frequency=f;
    allf=alltr(:,1:k);

    save([path 'spec__' name '.mat'], 'frequency','allf','averagef');
    disp('Terminated unregularly. Results so far saved.')
end