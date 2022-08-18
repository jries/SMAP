%vibration analysis playground
dir='/Users/ries/Library/Mobile Documents/com~apple~CloudDocs/mEMBL/Applications/Position/020_MaxPerutz/Negotiations/documents/vibrations/JonasRiesVibrations/Finale_Messung/';
files={};
files{1}={'MP1_0513.txt','MP2_0513.txt','MP3_0513.txt','MP4_0513.txt'};
files{2}={'MP5_4615.txt','MP6_4615.txt','MP7_4615.txt','MP8_4615.txt'};
files{2}={'MP9_4617.txt','MP10_4617.txt','MP11_4617.txt','MP12_4617.txt'};
files{3}={'MP13_2UG_8','MP14_2UG_8','MP15_2UG_8','MP16_2UG_8'};
files{4}={'MP17_2UG_8_2','MP18_2UG_8_2'};
files{5}={'MP19_2UG_8_1','MP20_2UG_8_1'};
files{6}={'MP21_2UG_8_3','MP22_2UG_8_3'};
% files{2}={'MP5_4615.txt','MP6_4615.txt','MP7_4615.txt','MP8_4615.txt'};
coord={'X','Y','Z'};
p.dtplot=5;

for c=1:length(coord)
    c
    p.linepar={'LineWidth',1};
    p.plottrace=true;
    figure(70+c);
    p.holdon=false;
    for k=1:length(files)
        plotvibrationanalysis(dir, files{k},coord{c},p);
        p.holdon=true;
    end
    p.linepar={'Color',[0 0 0],'LineWidth',2};
%     p.plottrace=false;
%     plotvibrationanalysis(dir, file,coord{c},p);
     legend({'0513','4615','4617','2UG_1','2UG_2','2UG_3','2UG_4'})
    subplot(2,3,2)
    title(coord{c})

    subplot(2,3,3)
    ax=gca;
    hold on
    loglog(ax.XLim,[1 1]*0.78)
    loglog(ax.XLim,[1 1]*1.56)
    loglog(ax.XLim,[1 1]*3.12)
    drawnow
    subplot(2,3,4)
    title(file)
end


%%
dir='/Users/ries/Data_local/ViennaVibrations/'; 
file='MP_0513_29_June.txt';
% file='MP_2UG_8_30Jun.txt';
file = 'MP_VBC2_2UG_19.4.txt';

% dir='/Volumes/Lacie/DataLacie/ViennaVibrations/'; 
% file='MP_VBC4_Raum_1105.txt';


dt=300; %s=5min
traw=readtable([dir, file]);
tx=1:dt:traw.Time(end);
ind1=1;
ttime={};
ttime{length(tx)-1}=traw(1:100,:);
for k=1:length(tx)-1
    while ind1<=length(traw.Time) && traw.Time(ind1)<tx(k)
        ind1=ind1+1;
    end
    ind2=ind1;
    while ind2<=length(traw.Time) && traw.Time(ind2)<tx(k+1)
        ind2=ind2+1;
    end
    ttime{k}=traw(ind1:ind2-1,:);
%     ttime{k}.Time=ttime{k}.Time-ttime{k}.Time(1);
    ind1=ind2;
end
files=ttime;
%%

coord={'X','Y','Z'};
p.dtplot=60;

for c=1:length(coord)
    c
    p.linepar={'LineWidth',1};
    p.plottrace=true;
    figure(89+c);
    p.holdon=false;
    for k=1:length(files)
        plotvibrationanalysis(dir, files{k},coord{c},p);
        p.holdon=true;
    end
    p.linepar={'Color',[0 0 0],'LineWidth',2};
    p.plottrace=false;
    plotvibrationanalysis(dir, file,coord{c},p);
%     legend({'0513','4615','2UG'})
    subplot(2,3,2)
    title(coord{c})

    subplot(2,3,3)
    ax=gca;
    hold on
    loglog(ax.XLim,[1 1]*0.78)
    loglog(ax.XLim,[1 1]*1.56)
    loglog(ax.XLim,[1 1]*3.12)
    drawnow
    subplot(2,3,4)
    title(file)
end

%%
files={};
files{1}='MP_0513_29_June.txt';
files{3}='MP_2UG_8_30Jun.txt';
files{2}='MP_VBC4_Raum_1105.txt';
files{4}='MP_VBC2_2UG_19.4.txt';
for c=1:length(coord)
    c
    p.linepar={'LineWidth',1};
    p.plottrace=true;
    figure(99+c);
    p.holdon=false;
    for k=1:length(files)
        plotvibrationanalysis(dir, files{k},coord{c},p);
        p.holdon=true;
    end
%     p.linepar={'Color',[0 0 0],'LineWidth',2};
%     p.plottrace=false;
%     plotvibrationanalysis(dir, file,coord{c},p);
%     legend({'VBC5_0513','VBC2_2UG','VBC4_1105'})
%     legend({'VBC5_0513','VBC4_1105','VBC2_2UG','VBC2_2UG_19'})
legend(files)
    subplot(2,3,2)
    title(coord{c})

    subplot(2,3,3)
    ax=gca;
    hold on
    loglog(ax.XLim,[1 1]*0.78)
    loglog(ax.XLim,[1 1]*1.56)
    loglog(ax.XLim,[1 1]*3.12)
    drawnow
%     subplot(2,3,4)
%     title(file)
end
% traw=readtable([dir file]);
% time=traw.Time;
% dt=0.001; %1 ms
% fdetect=1/dt;
% %play with Acc_Z;
% accz=traw.Acc_Z*1e6; %µm/s^2
% vz=cumsum(accz)*dt;
% % plot
% dtplot=0.1;
% tplot=0:dtplot:max(time);
% figure(88); subplot(2,2,1);
% plot(tplot,bindata(time,vz,tplot,'rms'))
% xlabel('time (s)')
% ylabel('v_{RMS} (µm/s)')
% 
% %psd
% subplot(2,2,2);
% hold off
% [vzPSD,fout]=periodogram(vz,[],[],fdetect);
% % loglog(fout,vzPSD,"Color",[1 1 1]*0.7)
% % hold on
% 
% df=fout(2)-fout(1);
% 
% % [azPSD,fout]=periodogram(accz,[],[],fdetect);
% % vzPSDa=azPSD./(2*pi*fout).^2;
% % loglog(fout,vzPSDa) %gives the same result!
% % [vzPSDs,fouts]=periodogram(vz(1:10000),[],[],1000);
% % loglog(fouts,vzPSDs)
% f_start=1;
% f_stop=512;
% octaves=log2(f_stop/f_start);
% 
% flog=logspace(log10(f_start),log10(f_stop),octaves*30+1);
% [vzPSDl,foutl]=periodogram(vz,[],flog,fdetect);
% subplot(2,2,2);
% 
% % loglog(foutl,vzPSDl,'k')
% 
% 
% %octave band
% f_to=logspace(log10(f_start),log10(f_stop),octaves*3+1);
% [vzm,fm]=bindata(fout,vzPSD,flog,'mean');
% loglog(flog,vzm)
% xlabel('frequency (Hz)')
% 
% [vztom,fwinm]=bindata(fout,vzPSD,f_to,'mean');
% % loglog(f_to,vztom,'b')
% 
% [vzto,foutto]=periodogram(vz,[],f_to,fdetect);
% subplot(2,2,2);
% ylabel('v PSD ((µm/s)^2/Hz)')
% % loglog(f_to,vzto,'c')
% 
% xlim([f_start,f_stop])
% ylim([10^-6 10^2])
% grid on
% 
% [vzbin,fwin]=bindata(fout,vzPSD,f_to,'sum');
% vrms=sqrt(vzbin.*df);
% subplot(2,2,3)
% hold off
% loglog(f_to(2:end-1), vrms(2:end-1))
% xlabel('frequency (Hz)')
% ylabel('v_{RMS} in 3rd octave (µm/s)')
% grid on
% % f_center=f_to(2:end-1)
% % f_edges=mean(vertcat(f_to(1:end-1),f_to(2:end)),1) %now on linear scale. I guess thats ok.
% % df_to=diff(f_edges)
% 
% subplot(2,2,4)
% semilogx(f_to(2:end-1), vrms(2:end-1))
% xlabel('frequency (Hz)')
% ylabel('v_{RMS} in 3rd octave (µm/s)')
% %v RMS octave band
% %1. use vzto and multiply by df, sqrt
% %2. use original and integrate over frequency band (as in Python)
% %  sum * df : we need df from before.
% %3. same also to calculate smooth curve.
% 
% %define frequency centers.
% %calculate edges
% %calculate df