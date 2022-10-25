%% clean up
clear all
close all

%% load and format localizations
SiR = readtable('SiR_localizations_for_fig4.csv');
mEos = readtable('mEos_localizations_for_fig4.csv');

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
load windows_for_fig4


%% generate reconstructed image
ref = default_iref([d(1).x(:),d(1).y(:)],25);
I1 = reconstruct(d(1), ref);
I2 = reconstruct(d(2), ref);
clear I
I(:, :, 1) = I1/1;
I(:, :, 2) = I2/4;
I(:, :, 3) = I1/1;
I = imgaussfilt(I, 2);

%% calculate the cross-correlations
dr = 50;
r = dr/2:dr:1000;
tau = frame_time:frame_time:2;

% basic cross-correlation, includes the edge correction;
[g,gerrs,time_edge_cor,N,Norm] = spacetime_xcor(d(1).x,d(1).y,d(1).t,d(2).x,d(2).y,d(2).t,spacewin,timewin,r,tau);
% gradient correction
[gm, gamma] = spatial_gradient_correction(d(1).x,d(1).y,d(2).x,d(2).y, spacewin, r, 1000);
g = g./gm;
gerrs = gerrs./gm;
% tau average for bleedthrough correction and improved statistics;
tauwindow = tau>.1  & tau <=2;
g_tauave = mean(g(:, tauwindow), 2);
dg_tauave = std(g(:, tauwindow),[], 2)/sqrt(sum(tauwindow));


%% make the figure
figure(1)
w1 = 6; h1=2.5;
set(1, 'Units', 'inches', 'color', 'w', 'PaperPosition', [0 0 w1 h1]);
pos = get(1, 'position');
set(1, 'position', [pos(1:2) w1 h1], 'Renderer', 'painters');

Dw = .08;
Dh = .15;
dw = .02;
dh = .02;
w = (1-3*Dw-dw)/3;
h = (1-1*Dh-dh)/1;

image_axis1 = [Dw-2*dw Dh-dh w+2*dw h+2*dh];
%image_axis2 = [Dw-2*dw Dh-2*dh w+2*dw h+2*dh];
xc_axis1 = [2*Dw+w Dh w h];
xc_axis2 = [3*Dw+2*w Dh w h];

alim = 1e4*[2.5 5.1 .7 3.4];

subplot('position', image_axis1)
imshow(I, ref)
hold on
plot(spacewin.p, 'facecolor', 'none', 'edgecolor', 'y')
plot(2.6e4+[0 5000], .9e4+[0 0], 'w-', 'linewidth', 3)
hold off
text(5e4, 3.1e4, {'BCR (SiR)'}, 'color', 'm', 'fontweight', 'bold', 'fontsize', 8, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right')
text(5e4, 3.3e4, {'GPI (mEos)'}, 'color', 'g', 'fontweight', 'bold', 'fontsize', 8, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right')

hold off
axis(alim)
axis off
text(2.2e4, .7e4, 'A', 'fontsize', 16)

%%

subplot('position', xc_axis1)
colors = colororder;
f = .5;
imax = 5;
alim = [0 2 .5 1.7];

% fit to gaussian plus offset for each tau
for i=1:imax
    subplot('position', [xc_axis1(1) sum(xc_axis1([2 4]))-xc_axis1(4)/imax*i  xc_axis1(3) xc_axis1(4)/imax])
    fill([tau tau(end:-1:1)], [g(i, :)+gerrs(i, :) g(i, end:-1:1)-gerrs(i, end:-1:1)], colors(i, :), 'edgecolor', 'none', 'facealpha', .2)
    hold on
    plot(tau, g(i, :),'.', 'color', colors(i, :));
    plot(tau, ones(size(tau)),'--', 'color',.5*[1 1 1]);
    
    F = fit(tau', g(i, :)', '1+A*exp(-x^2/2/s1^2)+B', 'startpoint', [.2 -.1 .1 ], 'lower', [0 -1 .05]);

    plot(tau, F(tau), 'color', colors(i, :));
    
    fitint(i) = F.B+1;
    ci = diff(confint(F, .68));
    dfitint(i) = ci(2)/2;
    
    axis(alim)
    text(1.9, .65, [num2str(r(i)) 'nm'], 'HorizontalAlignment', 'right', 'fontsize', 8)
    set(gca, 'ytick', .6:.4:1.4)
    if i<imax
        set(gca, 'xticklabel', [])
    else
        xlabel('separation time (sec)')
    end
    if i==3
        text(-.4, mean(alim(3:4)), 'cross-correlation', 'HorizontalAlignment', 'center', 'FontSize', 8, 'rotation', 90)
    end
    if i==1
        text(-.6, 1.7, 'B', 'fontsize', 16, 'VerticalAlignment', 'top')
    end
end

hold off

%%

subplot('position', xc_axis2)

H(1) = errorbar(r, g(:, 1), gerrs(:, 1), 'o-', 'markerfacecolor', colors(1, :));
hold on
H(3) = errorbar(r(1:imax), fitint,dfitint, 'v-', 'markerfacecolor', colors(2, :));
H(2) = errorbar(r, g_tauave, dg_tauave, 's-', 'markerfacecolor', colors(3, :), 'markersize', 3);
hold on
plot(r, ones(size(r)), '--', 'color', .5*[1 1 1])
hold off
axis([0 900 .9 1.5])


xlabel('separation distance (nm)')
ylabel('cross-correlation')
legend(H, 'g(\tau=0)', '<g(\tau<2s)>', 'g(\tau\rightarrow0)');
 text(-250, 1.5, 'C', 'fontsize', 16, 'VerticalAlignment', 'top')
print -dpng fig4.png