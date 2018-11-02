%% 
% GradientFit.m
% Used for 3D Emitter localization by Gradient direction Fitting
% This program is free software without any warranty. This set of programs 
% is distributed as part of academic resource sharing for non-profit 
% research only. Plaese see the GNU General Public License 
% for more details.
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Inputs
%   ROI  : Candidate emitter sub-region for position caculation
%          The center pixel defaults to the pixel with the peak intensity
%   RegR : Region radius of the ROI, so the size of ROI is (2*WinR)+1
%   GraR : The radius of Gradient used for caculation
%          Ignore the effect of the outer gradients with low SNR can
%          improve the performance of this method.        
% -------------------------------------------------------------------------
% Outputs
%   x,y  : The caculated emitter position in x-y dimensions
%          Unit: pixels
%   e    : The caculated ellipticity of the emitter's intensity distribution
% -------------------------------------------------------------------------
%
%
% Affiliation :  Departments of Medicine, University of Pittsburgh
% Reference   :  Ma, H. et al. Fast and Precise 3D Fluorophore Localization based on 
%                Gradient Fitting. Sci. Rep. 5, 14335; doi: 10.1038/srep14335 (2015).
% Link        :  http://www.pitt.edu/~liuy
%
% Copyright 2015, Hongqiang Ma
%%

function [x,y,e,wx,wy] = GradientFit3D(ROI,RegR,GraR)

% estimate the initial value of [x,y,e] by using centroid based method
[x0, y0, e0,wx,wy] = GradientFitCentroid(ROI);

% define the coordinates of the gradient gride, set the center pixel as the original point
[m,n] = meshgrid(0.5-GraR:GraR-0.5,GraR-0.5:-1:0.5-GraR);

% define the exact gradient at each position
Gx = e0*(x0-m);
Gy = (y0+n);
G2 = (Gx.^2 + Gy.^2);

% caculate the measured gradients
xID = RegR-GraR:RegR+GraR-1;
yID = xID;
gx = ROI(yID,xID+3)+ 2*ROI(yID+1,xID+3)+ 2*ROI(yID+2,xID+3)+ ROI(yID+3,xID+3)...
    -ROI(yID,xID)  - 2*ROI(yID+1,xID)  - 2*ROI(yID+2,xID)  - ROI(yID+3,xID);

gy = ROI(yID,xID)  + 2*ROI(yID,xID+1)  + 2*ROI(yID,xID+2)  + ROI(yID,xID+3)...
    -ROI(yID+3,xID)- 2*ROI(yID+3,xID+1)- 2*ROI(yID+3,xID+2)- ROI(yID+3,xID+3);

gx2 = gx.^2;
gy2 = gy.^2;
gxy = gx.*gy;
% g2 = (gx2 + gy2);

% Wg = sqrt(G2).*(g2);  % Wg is the weight parameter
% P = (G2.*g2)./Wg;
P = sqrt(G2);

% solve the equation to get the best fit [x,y,e] --------------------------
a1 = sum(sum(gy2.*m./P));
b1 = -sum(sum(gy2./P));
c1 = sum(sum(gxy./P));
d1 = -sum(sum(gxy.*n./P));

a2 = sum(sum(gxy.*m./P));
b2 = -c1;
c2 = sum(sum(gx2./P));
d2 = -sum(sum(gx2.*n./P));

a3 = sum(sum(gy2.*m.^2./P));
b3 = -2*a1;
c3 = -b1;
d3 = -d1;
e3 = -c1;
f3 = a2;
g3 = -sum(sum(gxy.*m.*n./P));

A1 = (a2*c1-a1*c2)/(b1*c2-b2*c1);
B1 = (c1*d2-c2*d1)/(b1*c2-b2*c1);

A2 = (a1*(b2*c1-b1*c2)+b1*(a1*c2-a2*c1))/(c1*(b1*c2-b2*c1));
B2 = (b1*(c2*d1-c1*d2)+d1*(b2*c1-b1*c2))/(c1*(b1*c2-b2*c1));

A = a3+A1*b3+A1*A1*c3+A1*A2*e3+A2*f3;
B = B1*B1*c3+B1*d3+B1*B2*e3;
C = B1*b3+2*A1*B1*c3+A1*d3+A1*B2*e3+B2*f3+g3;
%--------------------------------------------------------------------------

e = (-C+sqrt(C^2-4*A*B))/(2*A);
x = A1+B1/e;
y = A2*e+B2;
