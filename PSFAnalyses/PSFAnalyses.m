% Parameters
zstep = 20; % nm
pixsize = 130; % nm

% Import tiff
%[fileName,pathName] = uigetfile('*.tif')
srcFiles = dir(strcat(pathName,'\*.tif')); 
filename = strcat(pathName,srcFiles(1).name);
I(:,:,1) = imread(filename);
for i = 2 : length(srcFiles)
    filename = strcat(pathName,srcFiles(i).name);
    I(:,:,i) = imread(filename);
end

% Max intensity
maxI = squeeze(max(max(I)));
maxI_norm = maxI/max(maxI);

% Import table locs
loc = g.locData.loc;
frame = loc.frame;
x = loc.xnm;
y = loc.ynm;
PSFx = loc.PSFxnm;
PSFy = loc.PSFynm;

% Purge from frames with multiple localizations
for i=1:max(frame)
   ind = find(frame==i);
   if size(ind,1)>1
       frame(ind) = [];
       x(ind) = [];
       y(ind) = [];
       PSFx(ind) = [];
       PSFy(ind) = [];
   end
end

% Naive z0
ind_z0 = find(maxI == max(maxI));
z = [1:size(I,3)];
z = (z-ind_z0)*zstep; % nm and z0 = 0

% FWHM in z
FWHMz_maxI = FWHM(z,maxI);

% Find focus !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! need to fit with PSF
% ellipt not free
ind_focx = find(PSFx == min(PSFx));
ind_focy = find(PSFy == min(PSFy));

zfoc = z(frame(ind_focx));
focframe =frame(ind_focx);

% centroid x and y
meanx = mean(x);
meany = mean(y);
xc = x-meanx;
yc = y-meany;

xfoc = xc(ind_focx);
yfoc = yc(ind_focy);

r = sqrt((xc-xfoc).^2+(yc-yfoc).^2);

% FWHM x and y from image
width = size(I,1);
height = size(I,2);
frametot = size(I,3);
FWHMx = zeros(1,frametot);
FWHMy = zeros(1,frametot);
xprof = zeros(width,frametot);
yprof = zeros(height,frametot);
Ixprof = zeros(width,frametot);
Iyprof = zeros(height,frametot);

for k = 1:frametot
    maxrow = max(I(:,:,k));
    ind_maxyinfoc = find(maxrow == max(max(I(:,:,k))));
    maxcol = max(I(:,:,k)');
    ind_maxxinfoc = find(maxcol == max(max(I(:,:,k)')));
    
    if size(ind_maxyinfoc,2)>1 % find a better solution...
        ind_maxyinfoc = ind_maxyinfoc(1,1);
        ind_maxxinfoc = ind_maxxinfoc(1,1);
    end
    
    Iyprof(:,k) = I(ind_maxxinfoc,:,k);
    Ixprof(:,k) = I(:,ind_maxyinfoc,k);
    xprof(:,k) = ([1:size(Ixprof,1)]*pixsize)';
    yprof(:,k) = ([1:size(Iyprof,1)]*pixsize)';
    try
        FWHMx(k) = FWHM(xprof,Ixprof);
    catch
        FWHMx(k) = 0;
    end
    try
        FWHMy(k) = FWHM(yprof,Iyprof);
    catch
        FWHMy(k) = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Plots

% Max intensity with z
figure, plot(z,maxI)

% x,y scatter, z color-coded
figure, scatter(xc,yc,[],frame)

% distance to focus
figure, plot(frame,r)

% focus
figure, plot(z,FWHMx,'r.')
