function plotvibrationanalysis(dir, file,cor,p)
dt=0.001; %1 ms
if iscell(file)
    traw=readtable([dir, file{1}]);
    for k=2:length(file)
        th=readtable([dir, file{k}]);
        th.Time=th.Time+traw.Time(end)+dt;
        traw=vertcat(traw,th);
    end

elseif ischar(file)
traw=readtable([dir, file]);
elseif istable(file)
    traw=file;
end
if ~isfield(p,'linepar')
    p.linepar={};
end

time=traw.Time;

fdetect=1/dt;
%play with Acc_Z;
accz=traw.(['Acc_' cor])*1e6; %µm/s^2
vz=cumsum(accz)*dt;
% plot

for k=1:6
    subplot(2,3,k)
    if p.holdon
        hold on
    else
        hold off
    end
end



dtplot=p.dtplot;
tplot=0:dtplot:max(time);
 subplot(2,3,1);
 vzrms=bindata(time,vz,tplot,'std');
if p.plottrace
    plot(tplot,vzrms)
    % plot(tplot,bindata(time,vz,tplot,'rms'))
    xlabel('time (s)')
    ylabel('v_{std} (µm/s)')
else
    vzn=sqrt(mean(vzrms.^2));
    title(['Std(v): ' num2str(vzn,3) ' µm/s (' num2str(dtplot,3) ' s window)'])
end
%psd
subplot(2,3,2);

[vzPSD,fout]=periodogram(vz,[],[],fdetect);

df=fout(2)-fout(1);

f_start=1;
f_stop=512;
octaves=log2(f_stop/f_start);

flog=logspace(log10(f_start),log10(f_stop),octaves*30+1);
[vzPSDl,foutl]=periodogram(vz,[],flog,fdetect);


%octave band
f_to=logspace(log10(f_start),log10(f_stop),octaves*3+1);
[vzm,fm]=bindata(fout,vzPSD,flog,'mean');
loglog(flog,vzm,p.linepar{:})
xlabel('frequency (Hz)')

ylabel('v PSD ((µm/s)^2/Hz)')

xlim([f_start,f_stop])
ylim([10^-5 10^2])
grid on

[vzbin,fwin]=bindata(fout,vzPSD,f_to,'sum');
vrms=sqrt(vzbin.*df);
subplot(2,3,3)
% loglog(f_to([2 end-1]),[1 1]*0.78)
% hold on
% loglog(f_to([2 end-1]),[1 1]*1.56)
% loglog(f_to([2 end-1]),[1 1]*3.12)
% loglog(f_to([2 end-1]),[1 1]*6.25)
loglog(f_to(2:end-1), vrms(2:end-1),p.linepar{:})
xlabel('frequency (Hz)')
ylabel('v_{RMS} in 3rd octave (µm/s)')
grid on
% legend

subplot(2,3,4)
semilogx(f_to(2:end-1), vrms(2:end-1),p.linepar{:})
xlabel('frequency (Hz)')
ylabel('v_{RMS} in 3rd octave (µm/s)')

subplot(2,3,5);
if p.plottrace
    plot(tplot,bindata(time,accz,tplot,'std'))
end
% plot(tplot,bindata(time,vz,tplot,'rms'))
xlabel('time (s)')
ylabel('a_{std} (µm/s^2)')

subplot(2,3,6)
[azPSD,fout]=periodogram(accz,[],[],fdetect);
[azm,fm]=bindata(fout,azPSD,flog,'mean');
loglog(flog,azm,p.linepar{:})
xlabel('frequency (Hz)')

ylabel('acc PSD ((µm/s^2)^2/Hz)')
grid on
end