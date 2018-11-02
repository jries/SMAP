%% ****************************************************************************************************
%% Return 2D guassian kernel
%% 
%% Input:    x_inx,y_inx          : meshgrid of x and y 
%%           x_pos,y_pos          : center_position 
%%           Gsigma               : standard deviation of gaussian function
%%
%% Output:   img_kernel           : 2D Guassian kernel
%% ****************************************************************************************************

function img_kernel = Gauss_kernel(x_inx,y_inx,x_pos,y_pos,Gsigma)
    img_kernel = exp(-((x_inx-x_pos).^2+(y_inx-y_pos).^2)/(2*Gsigma^2))/(2*pi*Gsigma^2) ;
end