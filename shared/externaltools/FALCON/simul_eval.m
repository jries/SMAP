%% ****************************************************************************************************
%% calculation of localization error in x,y direction.
%% each localized particle is matched to the closest true particle within a radius
%% ****************************************************************************************************

function [num_idens,num_clusters,errors_x,errors_y] = simul_eval(Results,true_poses,pitch,num_frame,radius)
num_idens = zeros(num_frame,1);
num_clusters = zeros(num_frame,1);
errors_x = [];
errors_y = [];
for tt = 1:num_frame
    ind = find(Results(:,1) == tt);
    est_pos = Results(ind,2:3);
    true_pos = true_poses(:,:,tt);
    num_cluster = size(est_pos,1);
    if num_cluster
        num_iden_temp = zeros(1,size(true_pos,1));
        error_x = zeros(num_cluster,1);
        error_y = zeros(num_cluster,1);
        num_cluster_temp = zeros(num_cluster,1);
        for ii = 1: num_cluster
            % find the closest true location
            distance = sqrt((true_pos(:,1)-est_pos(ii,1)).^2 + (true_pos(:,2)-est_pos(ii,2)).^2);
            [error_2D,index] = min(distance*pitch);
            
            if abs(error_2D) < radius
                num_iden_temp(index) = 1;
                num_cluster_temp(ii) = 1;
                error_x(ii) = (true_pos(index,1)-est_pos(ii,1))*pitch;
                error_y(ii) = (true_pos(index,2)-est_pos(ii,2))*pitch;
            end
        end
        ind = find(num_cluster_temp>0);
        num_idens(tt) = sum(num_iden_temp(:));
        num_clusters(tt) = num_cluster;
        error_x = error_x(ind);
        error_y = error_y(ind);
        errors_x = [errors_x(:); error_x(:)];
        errors_y = [errors_y(:); error_y(:)];
    end
end
