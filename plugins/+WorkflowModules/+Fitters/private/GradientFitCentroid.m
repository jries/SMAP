% Estimate the initial value of [x,y,e] by using centroid based method
% -------------------------------------------------------------------------
% Inputs
%   ROI    : Candidate emitter sub-region for position estimation
% -------------------------------------------------------------------------
% Outputs
%   x0,y0  : The estimated emitter position in x-y dimensions
%          
%   e0     : The estimated ellipticity of the emitter's intensity distribution
% -------------------------------------------------------------------------

function [x0, y0, e0,Wx,Wy] = Centroid(ROI)
if sum(ROI(:))==0
    x0=0;y0=0; e0=0;Wx=0;Wy=0;
    return
end

[Ny, Nx] = size(ROI);
w = (Ny-1)/2;
X = 1:Nx;
Y = 1:Ny;
sIx = sum(ROI,1);
sIy = sum(ROI,2);
sIx = sIx - min(sIx(:));
sIy = sIy - min(sIy(:));
Xc = sum(sIx.*X)/sum(sIx);
Yc = sum(sIy.*Y')/sum(sIy);
X1 = 1-Xc:Nx-Xc;
Y1 = 1-Yc:Ny-Yc;
Wx = sum(sIx.*abs(X1))/sum(sIx);
Wy = sum(sIy.*abs(Y1'))/sum(sIy);

% for i=1:(2*w+1)
%     if i<Xc-3*Wx || i>Xc+3*Wx
%         sIx(i) = 0;
%     end
%     if i<Yc-3*Wy || i>Yc+3*Wy
%         sIy(i) = 0;
%     end
% end

sIx(1:floor(Xc-3*Wx)) = 0;
sIx(floor(Xc+3*Wx)+1:end) = 0;

sIx(1:floor(Yc-3*Wy)) = 0;
sIx(floor(Yc+3*Wy)+1:end) = 0;


Xc = sum(sIx.*X)/sum(sIx);
Yc = sum(sIy.*Y')/sum(sIy);
X1 = 1-Xc:Nx-Xc;
Y1 = 1-Yc:Ny-Yc;
Wx = sum(sIx.*abs(X1))/sum(sIx)/2; %two for correct radius, from radial center
Wy = sum(sIy.*abs(Y1'))/sum(sIy)/2;
x0 = Xc-w-1;
y0 = Yc-w-1;
e0 = Wy/Wx;

