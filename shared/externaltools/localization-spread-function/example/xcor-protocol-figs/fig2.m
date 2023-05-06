%% clean up
clear all
close all

%% load and format localizations
SiR = readtable('SiR_localizations_for_figs_2and3.csv');
mEos = readtable('mEos_localizations_for_figs_2and3.csv');

% only include times between 2 and 10 min after stimulation
keep1 = SiR.t>2*60 & SiR.t<10*60;
keep2 = mEos.t>2*60 & mEos.t<10*60;

% assemble data structure
d(1).x = SiR.x(keep1);
d(1).y = SiR.y(keep1);
d(1).t = SiR.t(keep1);
d(2).x = mEos.x(keep2);
d(2).y = mEos.y(keep2);
d(2).t = mEos.t(keep2);

frame_time = 0.0217; %in sec

% if needed, make spacewins (regions of interest) using spacewin_gui;
% also uses a timewin which includes breaks between sets of acquisitions.
load windows_for_fig2

%% generate reconstructed image
ref = default_iref([d(1).x(:),d(1).y(:)],25);
I1 = reconstruct(d(1), ref);
I2 = reconstruct(d(2), ref);
clear I
I(:, :, 1) = I1/.8;
I(:, :, 2) = I2/5;
I(:, :, 3) = I1/.8;
I = imgaussfilt(I, 2);

%% calculate the cross-correlations

dr = 50;
r = dr/2:dr:1000;
tau = frame_time:frame_time:2;

% internal mask
% basic cross-correlation, includes the edge correction;
[g,gerrs,time_edge_cor,N,Norm] = spacetime_xcor(d(1).x,d(1).y,d(1).t,d(2).x,d(2).y,d(2).t,spacewin_inside,timewin,r,tau);
% gradient correction
[gm, gamma] = spatial_gradient_correction(d(1).x,d(1).y,d(2).x,d(2).y, spacewin_inside, r, 1000);
g = g./gm;
gerrs = gerrs./gm;
% tau average for bleedthrough correction and improved statistics;
tauwindow = tau>.1  & tau <=2;
g_tauave_inside = mean(g(:, tauwindow), 2);
dg_tauave_inside = std(g(:, tauwindow),[], 2)/sqrt(sum(tauwindow));

% external mask
% basic cross-correlation, includes the edge correction;
[g,gerrs,time_edge_cor,N,Norm] = spacetime_xcor(d(1).x,d(1).y,d(1).t,d(2).x,d(2).y,d(2).t,spacewin_outside,timewin,r,tau);
% gradient correction
[gm, gamma] = spatial_gradient_correction(d(1).x,d(1).y,d(2).x,d(2).y, spacewin_outside, r, 1000);
g = g./gm;
gerrs = gerrs./gm;
% tau average for bleedthrough correction and improved statistics;
g_tauave_outside = mean(g(:, tauwindow), 2);
dg_tauave_outside = std(g(:, tauwindow),[], 2)/sqrt(sum(tauwindow));

%% make the figure
figure(1)
w1 = 4; h1=4;
set(1, 'Units', 'inches', 'color', 'w', 'PaperPosition', [0 0 w1 h1]);
pos = get(1, 'position');
set(1, 'position', [pos(1:2) w1 h1], 'Renderer', 'painters');

Dw = .1;
Dh = .1;
dw = .02;
dh = .02;
w = (1-2*Dw-dw)/2;
h = (1-2*Dh-dh)/2;

image_axis1 = [Dw-2*dw 2*Dh+h-2*dh w+2*dw h+2*dh];
image_axis2 = [Dw-2*dw Dh-2*dh w+2*dw h+2*dh];
xc_axis1 = [2*Dw+w 2*Dh+h w h];
xc_axis2 = [2*Dw+w Dh w h];

alim = 1e4*[2.5 5.1 .7 3.4];

barx = 2.6e4;
bary = .9e4;
textx = 2.6e4;
texty1 = 3.15e4;
texty2 = 3.35e4;

lettertexty = .7e4;
lettertextx = 2.1e4;

subplot('position', image_axis1)
imshow(I, ref)
hold on
plot(spacewin_outside.p, 'facecolor', 'none', 'edgecolor', 'y')
plot(barx+[0 5000], bary+[0 0], 'w-', 'linewidth', 3)
hold off
axis(alim)
axis off
text(textx, texty1, 'BCR (SiR)', 'color', 'm', 'fontweight', 'bold', 'fontsize', 8, 'VerticalAlignment', 'bottom')
text(textx, texty2, 'Src15 (mEos)', 'color', 'g', 'fontweight', 'bold', 'fontsize', 8, 'VerticalAlignment', 'bottom')
text(lettertextx,lettertexty, 'A', 'fontsize', 16,'VerticalAlignment', 'top')


subplot('position', image_axis2)
imshow(I, ref)
hold on
plot(spacewin_inside.p, 'facecolor', 'none', 'edgecolor', 'y')
plot(barx+[0 5000], bary+[0 0], 'w-', 'linewidth', 3)
hold off
axis(alim)
axis off
text(textx, texty1, 'BCR (SiR)', 'color', 'm', 'fontweight', 'bold', 'fontsize', 8, 'VerticalAlignment', 'bottom')
text(textx, texty2, 'Src15 (mEos)', 'color', 'g', 'fontweight', 'bold', 'fontsize', 8, 'VerticalAlignment', 'bottom')
text(lettertextx,lettertexty, 'B', 'fontsize', 16,'VerticalAlignment', 'top')
subplot('position', xc_axis1)


errorbar(r, g_tauave_outside, dg_tauave_outside, '.-')
hold on
plot(r, ones(size(r)), '--', 'color', .5*[1 1 1])
hold off
axis([0 900 .9 1.5])
xlabel('separation distance (nm)')
ylabel('cross-correlation')

subplot('position', xc_axis2)
errorbar(r, g_tauave_inside, dg_tauave_inside, '.-')
hold on
plot(r, ones(size(r)), '--', 'color', .5*[1 1 1])
hold off
axis([0 900 .9 1.1])
set(gca, 'ytick', [.9 1 1.1])

xlabel('separation distance (nm)')
ylabel('cross-correlation')

print -dpng fig2.png