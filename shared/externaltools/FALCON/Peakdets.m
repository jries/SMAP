%% find local peaks and get initial localizations by center of mass
function [est_c_new,delta_x,delta_y] = Peakdets(est_c,thresh)

[Ny,Nx,Nt] = size(est_c);
img = zeros(Ny,Nx,'single');
est_c_new = zeros(size(est_c),'single');
delta_x = zeros(size(est_c),'single');
delta_y = zeros(size(est_c),'single');

c_mask = ones(3);
x_mask = [-1 0 1;-1 0 1;-1 0 1];
x_mask = x_mask(:);
y_mask = [-1 -1 -1; 0 0 0; 1 1 1];
y_mask = y_mask(:);

est_c_copy = est_c;
est_c_copy(est_c<thresh) = 0;
for tt = 1: Nt
    img(2:end-1,2:end-1) = est_c_copy(2:end-1,2:end-1,tt);
    [I,J,~] = find(img>0);
    for ii = 1:length(I)
        xx = J(ii);
        yy = I(ii);
        img_temp = est_c(yy-1:yy+1,xx-1:xx+1,tt);
        
        % for local maxima
        if est_c(yy,xx,tt) >= max(img_temp(:))
            
            % center of mass
            img_temp2 = c_mask(:).*img_temp(:);
            photon = sum(img_temp2);
            dx = sum(x_mask.*img_temp2)/photon;
            dy = sum(y_mask.*img_temp2)/photon;
            pos_y = yy+round(dy);
            pos_x = xx+round(dx);
            est_c_new(pos_y,pos_x,tt) = photon;
            delta_x(pos_y,pos_x,tt) = dx-round(dx);
            delta_y(pos_y,pos_x,tt) = dy-round(dy);
        end
    end
end
end
