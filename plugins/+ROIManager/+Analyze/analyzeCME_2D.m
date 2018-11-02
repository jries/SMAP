global se;


r2D=0;
Ncirc=0;
si=0;
x0=0;
y0=0;
image_cumulative=0;

% for k=1:length(se.sites)
for k=1:153
        si=se.sites(k).evaluation.CME2DRingFit;
        r1(k)=si.r1;
        Ncirc(k)=si.Ncirc;
        x0(k)=si.x0;
        y0(k)=si.y0;
        image_cumulative=image_cumulative+se.sites(k).image.image;
end

% aligning the sites on center of fitted circle
if 0  
    for k=1:length(se.sites)
        
        se.sites(k).pos(1)=se.sites(k).pos(1)+x0(k);
        se.sites(k).pos(2)=se.sites(k).pos(2)+y0(k);
%         sitepar.sitelist.list(k).pos(1)=sitepar.sitelist.list(k).pos(1)+x0(k);
%         sitepar.sitelist.list(k).pos(2)=sitepar.sitelist.list(k).pos(2)+y0(k);
    end
end


% histogram of fitted 2D radii

figure(11)
clf
subplot(1,3,1)
histfit(r1,15)
[mu_r, sigma_r]=normfit(r1);
xlabel('radius');
title(strcat('Ede1: number of evaluated sites: ',sprintf('%d',length(r1))));

subplot(1,3,2)
histfit(Ncirc,15)
[mu_n, sigma_n]=normfit(Ncirc);
xlabel('number within 1.5x radius ??');
title(strcat('mu_r= ',sprintf('%.2f',mu_r),'+-',sprintf('%.2f',sigma_r),' mu_N= ',sprintf('%.2f',mu_n),'+-',sprintf('%.2f',sigma_n)));

subplot(1,3,3)
plot(Ncirc,r1,'.b')
xlabel('Ncirc'); ylabel('radius from fit');

% cumulative image



if 1
    figure(12)
    clf
    imagesc(mat2gray(imcrop(image_cumulative,[10 10 80 80])));
%     imagesc(mat2gray(image_cumulative));
    ax=gca;
%     ax.YLim = [0 400]; ax.XLim = [0 400];
    % image(se.sites(12).image.image);
end
