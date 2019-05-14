clc
close all
clearvars

% parameters
    fmin1 = 42;
    fmax1 = 52;

    fmin2 = 85;
    fmax2 = 105;

    fmin0 = fmin1;
    fmax0 = fmax2;

    gamma = 0.5;

    plot_frequencies = [43.7 48 48.9 49.4 50 98.1 100];
    plot_frequencies = sort(plot_frequencies);
% /parameters



[file path]=uigetfile([path filesep  '*.mat']);
%        file = 'spec__2019-05-13_1800__4pi_room__mic_pointing_at_Luxendo_wall__10s_1700blocks_ON.mat';
%        path = '/Users/diekmann/Documents/24 Audio spectra/';
load([path filesep file]);




% figure: spectra
figure
    for k = 1:length(plot_frequencies)
        % find index of corresponding frequency
            ind(k) = find( min( abs( frequency - plot_frequencies(k) ) ) == abs( frequency - plot_frequencies(k) ) );
    end % for k

    set(gcf, 'Position', [1 55 1440 750])

    % plot average spectra
    subplot(2,2,1)
        plot(frequency, averagef)
        hold on
        plot(frequency, prctile(allf',5))
        plot(frequency, max(allf'))
        legend('average', '5th percentile', 'max')
        xlim([fmin1 fmax1]); % fmin max 1
        xlabel('Frequency (Hz)')
        ylabel('Amp. (a.u.)')
        title(file)

    subplot(2,2,3)
        plot(frequency, averagef)
        hold on
        plot(frequency, prctile(allf',5))
        plot(frequency, max(allf'))
        legend('average', '5th percentile', 'max')
        xlim([fmin2 fmax2]); % fmin max 2
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

    
    
%figure: color representation
figure

    set(gcf, 'Position', [1 55 1440 750])

    subplot(1,2,1)
        ind_min0 = find( min( abs( frequency - fmin0 ) ) == abs( frequency - fmin0 ) );
        ind_max0 = find( min( abs( frequency - fmax0 ) ) == abs( frequency - fmax0 ) );

        imagesc(allf(ind_min0 : ind_max0, :).^gamma)
            xl = get(gca, 'XTick');
            yl = get(gca, 'YTick');
            set(gca, 'XTickLabel', cellstr(num2str(round(xl'*10/60))))
            set(gca, 'YTickLabel', cellstr(num2str(transpose(frequency(ind_min0+yl)))))

            xlabel('Time (min)')
            ylabel('Frequency (Hz)')
    
    % 3D plot of frequency range over time, window 1
    subplot(2,2,2)
        ind_min1 = find( min( abs( frequency - fmin1 ) ) == abs( frequency - fmin1 ) );
        ind_max1 = find( min( abs( frequency - fmax1 ) ) == abs( frequency - fmax1 ) );

        imagesc(allf(ind_min1 : ind_max1, :).^gamma)
            xl = get(gca, 'XTick');
            yl = get(gca, 'YTick');
            set(gca, 'XTickLabel', cellstr(num2str(round(xl'*10/60))))
            set(gca, 'YTickLabel', cellstr(num2str(transpose(frequency(ind_min1+yl)))))

            xlabel('Time (min)')
            ylabel('Frequency (Hz)')
            
    % 3D plot of frequency range over time, window 2
    subplot(2,2,4)
        ind_min2 = find( min( abs( frequency - fmin2 ) ) == abs( frequency - fmin2 ) );
        ind_max2 = find( min( abs( frequency - fmax2 ) ) == abs( frequency - fmax2 ) );

        imagesc(allf(ind_min2 : ind_max2, :).^gamma)
            xl = get(gca, 'XTick');
            yl = get(gca, 'YTick');
            set(gca, 'XTickLabel', cellstr(num2str(round(xl'*10/60))))
            set(gca, 'YTickLabel', cellstr(num2str(transpose(frequency(ind_min2+yl)))))

            xlabel('Time (min)')
            ylabel('Frequency (Hz)')