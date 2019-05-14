clc
%close all
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
set(gcf, 'Position', [1 55 1440 750])

% plot average spectrum
subplot(2,2,1)
    fmin = 85;
    fmax = 105;
    plot(frequency, averagef)
    hold on
    plot(frequency, prctile(allf',5))
    plot(frequency, max(allf'))
    legend('average', '5th percentile', 'max')
    xlim([fmin fmax]);
    xlabel('Frequency (Hz)')
    ylabel('Amp. (a.u.)')
    title(file)
   
% plot frequencies over time
for k = 1:length(plot_frequencies)
    subplot(length(plot_frequencies),2,2*k)
        plot((1:length(allf(ind(k), :)))*10/60, allf(ind(k), :))
        legend(num2str(frequency(ind(k))'))
        ylabel('Amp. (a.u.)')
end % for k
xlabel('time (minutes)')

% 3D plot of frequency range over time
subplot(2,2,3)
    ind_min = find( min( abs( frequency - fmin ) ) == abs( frequency - fmin ) );
    ind_max = find( min( abs( frequency - fmax ) ) == abs( frequency - fmax ) );

    imagesc(allf(ind_min : ind_max, :))
        xl = get(gca, 'XTick');
        yl = get(gca, 'YTick');
        set(gca, 'XTickLabel', cellstr(num2str(round(xl'*10/60))))
        set(gca, 'YTickLabel', cellstr(num2str(transpose(frequency(ind_min+yl)))))

        xlabel('Time (min)')
        ylabel('Frequency (Hz)')
    