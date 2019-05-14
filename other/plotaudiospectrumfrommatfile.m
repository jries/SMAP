clc
close all
clearvars

plot_frequencies = [43.7 48 48.9 49.4 50 98.1 100];
plot_frequencies = sort(plot_frequencies);

[file path]=uigetfile([path filesep  '*.mat']);
load([path filesep file]);


for k = 1:length(plot_frequencies)
    % find index of corresponding frequency
        ind(k) = find( min( abs( frequency - plot_frequencies(k) ) ) == abs( frequency - plot_frequencies(k) ) );
end % for k

figure
subplot(1,2,1)
    fmin = 30;
    fmax = 105;
    plot(frequency, averagef)
    hold on
    plot(frequency, min(allf'))
    plot(frequency, max(allf'))
    legend('average', 'min', 'max')
    xlim([fmin fmax]);
    xlabel('Frequency (Hz)')
    ylabel('Amplitude (a.u.)')
    title(file)

for k = 1:length(plot_frequencies)
    subplot(length(plot_frequencies),2,2*k)
        plot((1:length(allf(ind(k), :)))*10/60, allf(ind(k), :))
        legend(num2str(frequency(ind(k))'))
        ylabel('Amplitude (norm.)')
end % for k
xlabel('time (minutes)')

set(gcf, 'Position', [1 55 1440 750])