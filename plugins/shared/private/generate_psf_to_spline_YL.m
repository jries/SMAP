function [np_spline,coeff]=generate_psf_to_spline_YL(fname,s_size)
ndata = load(fname);
np_psf = double(ndata.psf);
% np_psf = np_psf/10000;
np_psf = np_psf/max(np_psf(:));
for i = 1:size(np_psf,3)
    np_psf1(:,:,i) = np_psf(:,:,i)/sum(sum(np_psf(:,:,13)));
end
% np_psf = np_psf/max(np_psf(:));
np_psf = np_psf1;
start = size(np_psf,1)/2-s_size;



s_size = 2*s_size;
np_spline = zeros(s_size,s_size,s_size);
xy_spline = [];
for i = 1:size(np_psf,3)
    xy_spline = [xy_spline Spline2D(np_psf(:,:,i))];
end

x = start;
for i = 1:s_size
    y = start;
    for j = 1:s_size
        zvals = zeros(size(np_psf,3),1);
%         zvals = zeros(51,1);
        for k = 1:size(np_psf,3)
            zvals(k) = xy_spline(k).f(y,x);
        end
        z_spline = Spline1D(zvals);

        max_z = size(np_psf,3)-1;
        inc = max_z/(s_size-1);
        for k = 1:s_size
            z = (k-1)*inc
            if (z> max_z)
                z = max_z;
            end
            np_spline(j,i,k) = z_spline.f(z);
        end
        y = y+1;
    end
    x = x+1;
end

% for i = 1:size(np_spline,3)
%     np_spline(:,:,i) = np_spline(:,:,i)/sum(sum(np_spline(:,:,i)));
% end

spline = Spline3D(np_spline);
coeff = spline.coeff;


