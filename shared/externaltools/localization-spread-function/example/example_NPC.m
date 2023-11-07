% Copyright (C) 2021 Thomas Shaw, and Sarah Veatch
% This file is part of SMLM SPACETIME RESOLUTION
% SMLM SPACETIME RESOLUTION is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% SMLM SPACETIME RESOLUTION is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% You should have received a copy of the GNU General Public License
% along with SMLM SPACETIME RESOLUTION.  If not, see <https://www.gnu.org/licenses/>

%% Resolution estimation example for a nuclear pore complex dataset

%load data
load('sampledata_npc.mat'); % The data is stored in a struct called 'data'
data

%% Optionally draw a new spatial window / ROI. Skip this step to use the one from the paper.
% spacewin_gui is a helper gui for drawing spatial windows that may have
% holes or multiple disjoint segments. Press 'save and close' when
% you are done.
% data.spacewin = spacewin_gui(data, 'PixelSize', 10) % use 10nm pixels

%% Run the resolution estimation routine
% Here we supply the data as a struct with fields x,y,t,spacewin,timewin.
% The function also accepts these arguments separately (in that order).
[corrdata, params] = spacetime_resolution(data, 'NTauBin', 10, 'Bootstrap', true);

%% Plot the correlation functions for each tau bin
figure;
plot(corrdata.r, corrdata.cWA)
tau = corrdata.taubincenters;
lh = legend(arrayfun(@num2str, tau, 'UniformOutput', false));
title(lh,'\tau (s)');
xlabel 'r (nm)'
ylabel 'g(r, \tau)'
set(gca, 'YScale', 'log');

%% Plot the normalized correlation function differences
figure;
plot(corrdata.r, corrdata.nDg);
lh = legend(arrayfun(@num2str, tau(1:end-1), 'UniformOutput', false));
title(lh,'\tau (s)');
xlabel 'r (nm)'
ylabel 'g(r, \tau)'

%% Plot the estimated resolution as a function of tau, and report average resolution
figure;
errorbar(tau(1:end-1),corrdata.s,corrdata.confint, 'o-');
title(sprintf('Average resolution is %.1f nm ', corrdata.S))
xlabel '\tau (s)'
ylabel 'resolution estimate (nm)'
