%% rearrange particles having high displacement values
function [est_c_new,delta_x_new,delta_y_new,delta_z_new] = handover(est_c,delta_x,delta_y,delta_z)
% find indices of particles which are close to adjacent grid points.
ind_arrange = find((abs(delta_x)>0.99)|(abs(delta_y)>0.99));
est_c_new = est_c;
delta_x_new = delta_x;
delta_y_new = delta_y;
delta_z_new = delta_z;

% move the particles to the closest grid points
% and adjust displacement terms.

if ~isempty(ind_arrange)
    [Ny,Nx,Nz] = size(est_c);
    est_c_new(ind_arrange)=0;
    delta_x_new(ind_arrange)=0;
    delta_y_new(ind_arrange)=0;
    delta_z_new(ind_arrange)=0;
    
    weight = est_c_new>0;
    [pos_y,pos_x,pos_z] = ind2sub(size(est_c),ind_arrange);
    pos_y_new = min(max(pos_y + delta_y(ind_arrange),0.51),Ny+0.49);
    pos_x_new = min(max(pos_x + delta_x(ind_arrange),0.51),Nx+0.49);
    
    ind_new = sub2ind([Ny,Nx,Nz],round(pos_y_new),round(pos_x_new),pos_z);
    
    [ind_new2,ia] = unique(gather(ind_new));
    delta_x_new_temp = pos_x_new - round(pos_x_new);
    delta_y_new_temp = pos_y_new - round(pos_y_new);
    delta_z_new_temp = delta_z(ind_arrange);
    photon_temp = est_c(ind_arrange);
    
    delta_x_new(ind_new2) = delta_x_new(ind_new2) + delta_x_new_temp(ia);
    delta_y_new(ind_new2) = delta_y_new(ind_new2) + delta_y_new_temp(ia);
    delta_z_new(ind_new2) = delta_z_new(ind_new2) + delta_z_new_temp(ia);
    est_c_new(ind_new2) = est_c_new(ind_new2)+photon_temp(ia);
    weight(ind_new2) = weight(ind_new2)+1;
    
    % In case that two particles are merged into the same grid bin
    if length(ind_new2)<length(ind_new)
        ind_new(ia) = 0;
        ind_temp = ind_new>0;
        delta_x_new_temp = delta_x_new_temp(ind_temp);
        delta_y_new_temp = delta_y_new_temp(ind_temp);
        delta_z_new_temp = delta_z_new_temp(ind_temp);
        photon_temp = photon_temp(ind_temp);
        [ind_new3,ia] = unique(gather(ind_new(ind_temp)));
        delta_x_new(ind_new3) = delta_x_new(ind_new3) + delta_x_new_temp(ia);
        delta_y_new(ind_new3) = delta_y_new(ind_new3) + delta_y_new_temp(ia);
        delta_z_new(ind_new3) = delta_z_new(ind_new3) + delta_z_new_temp(ia);
        est_c_new(ind_new3) = est_c_new(ind_new3)+photon_temp(ia);
        weight(ind_new3) = weight(ind_new3)+1;
    end
    
    % normalization
    delta_x_new(ind_new2) = delta_x_new(ind_new2)./weight(ind_new2);
    delta_y_new(ind_new2) = delta_y_new(ind_new2)./weight(ind_new2);
    delta_z_new(ind_new2) = delta_z_new(ind_new2)./weight(ind_new2);
end


end

