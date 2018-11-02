function [Results] = Peakdet2D(img,thresh,radius)
Results = [];
c_mask = ones(3);
x_mask = [-1 0 1;-1 0 1;-1 0 1];
x_mask = x_mask(:);
y_mask = [-1 -1 -1; 0 0 0; 1 1 1];
y_mask = y_mask(:);

img_copy = zeros(size(img),'single');
img_copy((radius+1):end-radius,(radius+1):end-radius) =img((radius+1):end-radius,(radius+1):end-radius);
img_copy(img_copy<thresh) = 0;
[I,J] = find(img_copy>0);
Results = zeros(length(I),3,'single');
for ii = 1:length(I)
    xx = J(ii);
    yy = I(ii);
    img_temp = img(yy-radius:yy+radius,xx-radius:xx+radius);
    img_temp_boundary = img_temp;
    img_temp_boundary(2:end-1,2:end-1) = 0;
    % for local maxima
    %     if (img(yy,xx) == max(img_temp(:)))&& (img(yy,xx)>1.2*max(img_temp_boundary(:)))
    
    if img(yy,xx) == max(img_temp(:))
        img_temps = img_temp(radius:(radius+2),radius:(radius+2));
        % center of mass
        img_temp2 = c_mask(:).*img_temps(:);
        photon_temp = sum(img_temp2);
        dx = sum(x_mask.*img_temp2)/photon_temp;
        dy = sum(y_mask.*img_temp2)/photon_temp;
        pos_y_temp = yy+round(dy);
        pos_x_temp = xx+round(dx);
        Results_temp = [pos_x_temp,pos_y_temp,photon_temp];
        Results = [Results;Results_temp];
    end
end
ind = Results(:,3)>0;
Results = Results(ind,:);
end
