% function [delta_f,delta_dxf,delta_ddxf,delta_dyf,delta_ddyf,delta_dzf,delta_ddzf] = computeDelta3Dj(z_delta, y_delta, x_delta)
function [delta_f,delta_dxf,delta_dyf,delta_dzf] = computeDelta3Dj_vec(x_delta, y_delta, z_delta)

% if x_delta<0
%     x_delta = 0;
% end
% 
% if x_delta>1
%     x_delta = 1;
% end
% 
% if y_delta<0
%     y_delta = 0;
% end
% 
% if y_delta>1
%     y_delta = 1;
% end
% 
% if z_delta<0
%     z_delta = 0;
% end
% 
% if z_delta>1
%     z_delta = 1;
% end

delta_f = (zeros(1,1,1,64,1)); %64 is number of basis functions (for example, xyz or x^2yz)
delta_dxf = (zeros(1,1,1,64,1));
delta_dyf = (zeros(1,1,1,64,1));
delta_dzf = (zeros(1,1,1,64,1));
cz = (1);
for i = 0:3
    cy = (1);
    for j = 0:3
        cx = (1);
        for k = 0:3
            delta_f(i*16+j*4+k+1) = (cx*cy*cz);
            if k<3
                delta_dxf(i*16+j*4+k+1+1) = ((k+1))*cx*cy*cz;
            end

            if j<3
                delta_dyf(i*16+(j+1)*4+k+1) = ((j+1))*cx*cy*cz;
            end

            if i<3
                delta_dzf((i+1)*16+j*4+k+1) = ((i+1))*cx*cy*cz;
            end

            cx = cx*x_delta;
        end
        cy = cy*y_delta;
    end
    cz = cz*z_delta;
end