%% clean up
clear all
close all

%% load and format localizations
SiR = readtable('SiR_localizations_for_figs_2and3.csv');
mEos = readtable('mEos_localizations_for_figs_2and3.csv');

% only include times before stimulation
keep1 = SiR.t<0;
keep2 = mEos.t<0;

% assemble data structure
d(1).x = SiR.x(keep1);
d(1).y = SiR.y(keep1);
d(1).t = SiR.t(keep1);
d(2).x = mEos.x(keep2);
d(2).y = mEos.y(keep2);
d(2).t = mEos.t(keep2);

frame_time = 0.0217; %in sec

% load spacewins (regions of interest), which were made using spacewin_gui;
% and a timewin which includes breaks between sets of acquisitions.
load windows_for_fig2

% If you'd like to try a different spacewin, run the following:
% spacewin_inside = spacewin_gui(d);

%% make reconstructed images for display
ref = default_iref([d(1).x(:),d(1).y(:)],25);
I1 = reconstruct(d(1), ref);
I2 = reconstruct(d(2), ref);

clear I
I(:, :, 1) = I1/1.5;
I(:, :, 2) = I2/2;
I(:, :, 3) = I1/1.5;
I = imgaussfilt(I, 2);

I1rgb = zeros(size(I));
I1rgb(:, :, 1) = I1/1.5;
I1rgb(:, :, 3) = I1/1.5;
I1rgb = imgaussfilt(I1rgb, 2);

I2rgb = zeros(size(I));
I2rgb(:, :, 2) = I2/2;
I2rgb = imgaussfilt(I2rgb, 2);

%% calculate cross-correlations

dr = 50;
r = dr/2:dr:1000; % centers of R bins in nm
tau = frame_time:frame_time:2; % centers of tau bins in sec

% basic cross-correlation, includes the edge correction, but no density gradient correction;
[g,gerrs,time_edge_cor,N,Norm] = spacetime_xcor(d(1).x,d(1).y,d(1).t,...% first data set
    d(2).x,d(2).y,d(2).t,...% second data set
    spacewin_inside,timewin,...%spatial and temporal windows/region of interest
    r,tau); % distance and time separations to compute cross-correlation at

% gradient correction
[gm, gamma] = spatial_gradient_correction(d(1).x,d(1).y,d(2).x,d(2).y, spacewin_inside, r, 1000);
gc = g./gm;
gcerrs = gerrs./gm;

% tau average for bleedthrough correction and improved statistics;
% do for both gradient corrected and uncorrected versions
tauwindow = tau>.1  & tau <=2; % in sec, i.e. consider localizations between 100 ms and 2s apart in time
g_tauave = mean(g(:, tauwindow), 2);
dg_tauave = std(g(:, tauwindow),[], 2)/sqrt(sum(tauwindow));
gc_tauave = mean(gc(:, tauwindow), 2);
dgc_tauave = std(gc(:, tauwindow),[], 2)/sqrt(sum(tauwindow));

% this calculates the edge correction explicitly.  Called inside
% spacetime_xcor but needed for this figure.
edge_cor = spatial_edge_correction(spacewin_inside, r);

%% make the figure

figure(1)
w1 = 4; h1=3;
set(1, 'Units', 'inches', 'color', 'w', 'PaperPosition', [0 0 w1 h1]);
pos = get(1, 'position');
set(1, 'position', [pos(1:2) w1 h1], 'Renderer', 'painters');

Dw = .1;
Dh = .15;
dw = .02;
dh = .02;
w = (1-Dw-3*dw)/3;
h = (1-Dh-2*dh)/2;

image_axis1 = [Dw Dh+dh+h w h];
image_axis2 = [Dw+w+dw Dh+dh+h w h];
image_axis3 = [Dw+2*w+2*dw Dh+dh+h w h];

xc_axis1 = [Dw Dh w h];
xc_axis2 = [Dw+w+dw Dh w h];
xc_axis3 = [Dw+2*w+2*dw Dh w h];

alim = 1e4*[2.5 5.1 .7 3.4];
barx = 2.6e4;
bary = .9e4;
textx = 2.6e4;
texty = 3.3e4;
lettertexty = .7e4;
lettertextx = 1.7e4;
lettertexty2 = 3.4e4;

subplot('position', image_axis1)
imshow(I1rgb, ref)
hold on
plot(spacewin_inside.p, 'facecolor', 'none', 'edgecolor', 'y')
plot(barx+[0 5000], bary+[0 0], 'w-', 'linewidth', 3)
hold off
text(textx, texty, {'BCR' '(SiR)'}, 'color', 'w', 'fontweight', 'bold', 'fontsize', 8, 'VerticalAlignment', 'bottom')
axis(alim)
axis off
text(lettertextx,lettertexty, 'A', 'fontsize', 16,'VerticalAlignment', 'top')
text(lettertextx,lettertexty2, 'B', 'fontsize', 16,'VerticalAlignment', 'top')

subplot('position', image_axis2)
imshow(I2rgb, ref)
hold on
plot(spacewin_inside.p, 'facecolor', 'none', 'edgecolor', 'y')
plot(barx+[0 5000], bary+[0 0], 'w-', 'linewidth', 3)
hold off
text(textx, texty, {'Src15' '(mEos3.2)'}, 'color', 'w', 'fontweight', 'bold', 'fontsize', 8, 'VerticalAlignment', 'bottom')
axis(alim)
axis off

subplot('position', image_axis3)
imshow(I, ref)
hold on
plot(spacewin_inside.p, 'facecolor', 'none', 'edgecolor', 'y')
plot(barx+[0 5000], bary+[0 0], 'w-', 'linewidth', 3)
hold off
text(textx, texty, {'merge'}, 'color', 'w', 'fontweight', 'bold', 'fontsize', 8, 'VerticalAlignment', 'bottom')

axis(alim)
axis off

subplot('position', xc_axis1)

textx = 50;
texty = .855;
alim = [0 900 .85 1.1];

errorbar(r, g_tauave.*edge_cor, dg_tauave.*edge_cor, '.-')
hold on
plot(r, ones(size(r)), '--', 'color', .5*[1 1 1])
hold off
axis(alim)
xlabel('distance (nm)')
ylabel('c(r)')
text(textx, texty, 'no corrections', 'HorizontalAlignment', 'left', 'fontsize', 8, 'VerticalAlignment', 'bottom')


subplot('position', xc_axis2)
errorbar(r, g_tauave, dg_tauave, '.-')
hold on
plot(r, ones(size(r)), '--', 'color', .5*[1 1 1])
hold off
axis(alim)
set(gca, 'yticklabel', [])
text(textx, texty, 'edge correction', 'HorizontalAlignment', 'left', 'fontsize', 8, 'VerticalAlignment', 'bottom')
xlabel('distance (nm)')

subplot('position', xc_axis3)
errorbar(r, gc_tauave, dgc_tauave, '.-')
hold on
plot(r, ones(size(r)), '--', 'color', .5*[1 1 1])
hold off
axis(alim)
set(gca, 'yticklabel', [])

xlabel('distance (nm)')
text(textx, texty, {'edge and' 'gradient' 'corrections'}, 'HorizontalAlignment', 'left', 'fontsize', 8, 'VerticalAlignment', 'bottom')