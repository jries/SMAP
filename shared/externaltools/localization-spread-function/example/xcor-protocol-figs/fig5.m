%% clean up
clear all
close all

%% load and format localizations
SiR = readtable('SiR_localizations_for_fig4.csv');
mEos = readtable('mEos_localizations_for_fig4.csv');

% only include times before stimulation
keep = mEos.t>2*60 & mEos.t<10*60;

% assemble data structure
d(1).x = mEos.x(keep);
d(1).y = mEos.y(keep);
d(1).t = mEos.t(keep);

ind = find(d(1).t>d(1).t(1), 1);
frame_time = d(1).t(ind)-d(1).t(1);
integration_time = .02; %sec

% if needed, make spacewins (regions of interest) using spacewin_gui;
% also uses a timewin which includes breaks between sets of acquisitions.
load windows_for_fig4

ref = default_iref([d(1).x(:),d(1).y(:)],25);
I1 = reconstruct(d(1), ref);
clear I
I(:, :, 2) = I1/3;
I(:, :, 1) = 0;
I(:, :, 3) = 0;
I = imgaussfilt(I, 2);

%% calculate the auto-correlation
dr = 20;
r = dr/2:dr:1000;
tau = frame_time:frame_time:2;

% basic cross-correlation, includes the edge correction;
[g,gerrs,time_edge_cor,N,Norm] = spacetime_xcor(d(1).x,d(1).y,d(1).t,d(1).x,d(1).y,d(1).t,spacewin,timewin,r,tau);

% gradient correction
[gm, gamma] = spatial_gradient_correction(d(1).x,d(1).y,d(1).x,d(1).y, spacewin, r, 1000);
g = g./gm;
gerrs = gerrs./gm;

%% fit to extract MSDs and population ratios
g_one_MSD = fittype('1/4/s1^2*N*1e6*exp(-x.^2/4/s1^2)');
g_two_MSD = fittype('alpha/4/s1^2*N*1e6*exp(-x.^2/4/s1^2)+(1-alpha)/4/s2^2*N*1e6*exp(-x.^2/4/s2^2)');

taumax = .75;
startpoint = [.1 1 30 100];
for i=1:numel(tau(tau<taumax))
    [rr, gg] = prepareCurveData(r', g(:, i));
    F_g = fit(rr, gg-1, g_two_MSD, 'startpoint', startpoint);
    MSD_slow(i) = 4*F_g.s1^2;
    MSD_fast(i) = 4*F_g.s2^2;
    alpha(i) = F_g.alpha;
    F_gs{i} = F_g;
end

tau_corrected = tau(tau<taumax).*(1-1/3*integration_time./tau(tau<taumax));

%% make the figure
figure(1)
w1 = 5; h1=5;
set(1, 'Units', 'inches', 'color', 'w', 'PaperPosition', [0 0 w1 h1]);
pos = get(1, 'position');
set(1, 'position', [pos(1:2) w1 h1], 'Renderer', 'painters');

Dw = .15;
Dh = .1;
dw = .04;
dh = .02;
w = (1-2*Dw-dw)/2;
h = (1-2*Dh-dh)/2;
h2 =(h)/3;
h3 =(h)/2;

image_axis1 = [Dw-2*dw 2*Dh+h w+2*dw h];
ac_axis = [2*Dw+w 2*Dh+h w h];
MSD_fast_axis = [Dw Dh+2*h2 w h2];
MSD_slow_axis = [Dw Dh+h2 w h2];
alpha_axis = [Dw Dh w h2];

D_axis = [2*Dw+w Dh+h3 w h3];
CR_axis = [2*Dw+w Dh w h3];

alim = 1e4*[2.5 5.1 .7 3.4];

subplot('position', image_axis1)
imshow(I, ref)
hold on
plot(spacewin.p, 'facecolor', 'none', 'edgecolor', 'y')
plot(2.6e4+[0 5000], .9e4+[0 0], 'w-', 'linewidth', 3)
hold off
text(5e4, 3.3e4, {'GPI (mEos)'}, 'color', 'w', 'fontweight', 'bold', 'fontsize', 8, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right')

hold off
axis(alim)
axis off
text(2.1e4, .7e4, 'A', 'fontsize', 16, 'VerticalAlignment', 'top')

subplot('position', ac_axis)
colors = colororder;
f = .5;
imax = 5;
xlim = [0 1];

% plot fits at the first few taus
for i=1:imax
    subplot('position', [ac_axis(1) sum(ac_axis([2 4]))-ac_axis(4)/imax*i  ac_axis(3) ac_axis(4)/imax])
    fill([r r(end:-1:1)], [g(:, i)'+gerrs(:, i)' g(end:-1:1, i)'-gerrs(end:-1:1, i)'], colors(i, :), 'edgecolor', 'none', 'facealpha', .2)
    hold on
    plot(r, g(:, i),'.', 'color', colors(i, :));
    plot(r, ones(size(r)),'--', 'color',.5*[1 1 1]);
    
    plot(r, F_gs{i}(r)+1, 'color', colors(i, :));
    
    ylim([-.1*g(1, i) 1.3*g(1, i)])
    text(900, 1*g(1, i), [num2str(tau_corrected(i), 2) 's'], 'HorizontalAlignment', 'right', 'fontsize', 8)
    if g(1, i) >50
        set(gca, 'ytick', 0:50:1*g(1, i))
    elseif g(1, i) >20
        set(gca, 'ytick', 0:20:1*g(1, i))
    else
        set(gca, 'ytick', 0:10:1*g(1, i))
    end
    
    if i<imax
        set(gca, 'xticklabel', [])
    else
        xlabel('separation distance (nm)')
    end
    if i==3
        text(-200, mean(get(gca, 'ylim')), 'auto-correlation', 'HorizontalAlignment', 'center', 'FontSize', 8, 'rotation', 90)
    end
    if i==1
        text(-250, max(get(gca, 'ylim')), 'B', 'fontsize', 16, 'VerticalAlignment', 'top')
    end
end

hold off

subplot('position', MSD_fast_axis)
plot(tau_corrected, MSD_fast*1e-6, 'o', 'markersize', 3)
ylim([0 1])
text(-.1, mean(get(gca, 'ylim')), {'fast MSD' '(\mum^2)'}, ...
    'rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'fontsize', 9)

set(gca, 'xticklabel', [])
hold off
text(-.3, 1.3, 'C', 'fontsize', 16, 'VerticalAlignment', 'top')


subplot('position', MSD_slow_axis)
plot(tau_corrected, MSD_slow/1000, 'o', 'markersize', 3)
set(gca, 'xticklabel', [])
ylim([1.4 4.9])
text(-.1, mean(get(gca, 'ylim')), {'slow MSD' '(10^3nm^2)'}, ...
    'rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'fontsize', 9)

subplot('position', alpha_axis)
plot(tau_corrected, alpha*100, 'o', 'markersize', 3)
ylim([1 3.5])
xlabel('separation time (sec)')
text(-.1, mean(get(gca, 'ylim')), {'% slow'}, ...
    'rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'fontsize', 9)


subplot('position', D_axis)
LP = 30;
plot(tau_corrected, (MSD_fast-2*LP^2)*1e-6/4./tau_corrected, 'o', 'markersize', 3)
ylim([.15 .5])
set(gca, 'xticklabels', [])
text(-.1, mean(get(gca, 'ylim')), {'diffusion coeff.' '(\mum^2/sec)'}, ...
    'rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'fontsize', 9)
text(-.27, .58, 'D', 'fontsize', 16, 'VerticalAlignment', 'top')

subplot('position', CR_axis)
plot(tau_corrected, sqrt((MSD_slow-2*LP^2)), 'o', 'markersize', 3)
ylim([0 60])
xlabel('separation time (sec)')
text(-.1, mean(get(gca, 'ylim')), {'conf. radius' '(nm)'}, ...
    'rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'fontsize', 9)