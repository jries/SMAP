classdef CME2CSide_yule2<interfaces.SEEvaluationProcessor
    properties
        savedevals
    end
    methods
        function obj=CME2CSide_yule2(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function makeGui(obj)
            makeGui@interfaces.SEEvaluationProcessor(obj);
%             obj.guihandles.saveimagesb.Callback={@saveimagesb_callback,obj};
        end
        function out=run(obj,p)
            out=runintern(obj,p);
         
        end
        
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end

function pard=guidef

pard.t1.object=struct('Style','text','String','segment:');
pard.t1.position=[1,1];
pard.t1.Width=1;
pard.segmentl1.object=struct('Style','checkbox','String','layer1','Value',1);
pard.segmentl1.position=[1,2];
pard.segmentl2.object=struct('Style','checkbox','String','layer2','Value',1);
pard.segmentl2.position=[1,3];
pard.segmentl1l2.object=struct('Style','checkbox','String','l1+l2','Value',1);
pard.segmentl1l2.position=[1,4];

pard.t2.object=struct('Style','text','String','layer1');
pard.t2.position=[4,3];
pard.t2.Width=1;
pard.t3.object=struct('Style','text','String','layer2');
pard.t3.position=[4,4];
pard.t3.Width=1;
pard.t4.object=struct('Style','text','String','filter [sx sz]');
pard.t4.position=[5,1];
pard.t4.Width=2;
pard.sigma1.object=struct('Style','edit','String','15 5');
pard.sigma1.position=[5,3];
pard.sigma2.object=struct('Style','edit','String','5 5');
pard.sigma2.position=[5,4];

pard.t5.object=struct('Style','text','String','cutoff factor');
pard.t5.position=[2,1];
pard.t5.Width=2;
pard.cutofffactor.object=struct('Style','edit','String','1');
pard.cutofffactor.position=[2,3.5];

pard.t6.object=struct('Style','text','String','dilation');
pard.t6.position=[3,1];
pard.t6.Width=2;
pard.dilation.object=struct('Style','edit','String','2');
pard.dilation.position=[3,3.5];

pard.plugininfo.type='ROI_Evaluate';

pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi','layer1_','layer2_','se_sitepixelsize'};
end



function out=runintern(obj,p)
if length(p.sigma1)==1
    p.sigma1(2)=p.sigma1(1);
end
if length(p.sigma2)==1
    p.sigma2(2)=p.sigma2(1);
end
roisize=ones(2,1)*p.se_siteroi/2;
roisizeh=p.se_siteroi/2;

[mask,mask1,mask2,mask12,imf1,imf2]=makemaskall(p,obj);
locs1=obj.getLocs({'xnm','ynm','locprecnm','xnmrot','ynmrot'},'layer',1,'size',roisizeh);
locs2=obj.getLocs({'xnm','ynm','locprecnm','xnmrot','ynmrot'},'layer',2,'size',roisizeh);

out.locs1 = locs1;
out.locs2 = locs2;


% use contour to segment
if 0
    pos1 = [];
    pos1.x = locs1.xnmrot;
    pos1.y = locs1.ynmrot;
    pos1.sx = zeros(size(pos1.x),'like',pos1.x);
    pos1.sy = zeros(size(pos1.x),'like',pos1.y);
    pos1.sx = pos1.sx + 12; 
    pos1.sy = pos1.sy + 12; 
    pos1.c = 0;


    premask1 = gaussrender_ellipt(pos1, [-roisizeh roisizeh],  [-roisizeh roisizeh], 1, 1,0, [0 0],25);
    [~,h1] = contour(premask1,100);
    h1Contour75 = h1.LevelList(15);
    [h1C75,~] = contour(premask1,[h1Contour75 h1Contour75]);
    h1C75x = floor(h1C75(1,2:(end-1)));
    h1C75y = floor(h1C75(2,2:(end-1)));
    theMask = zeros(300);
    theMask(sub2ind([300 300], h1C75y+1, h1C75x+1)) = 1;
    theMask = imfill(theMask,'holes');
    figure; imagesc(theMask)
    %figure; contour(premask1,100);
    hold on
    scatter(pos1.x+150, pos1.y+150)
end

% locs1.xnmrot = PCArot1(:,1); locs1.ynmrot = PCArot1(:,2);
% locs2.xnmrot = PCArot2(:,1); locs2.ynmrot = PCArot2(:,2);



inmask1=inmask(p,locs1,mask);
inmask2=inmask(p,locs2,mask);

xm1=locs1.xnmrot(inmask1);
ym1=locs1.ynmrot(inmask1);
xm2=locs2.xnmrot(inmask2);
ym2=locs2.ynmrot(inmask2);
%determine parameters

%% do PCA on the coordinates
if 0
    % Do PCA together
    oriCor = [xm1 ym1;xm2 ym2];
    [U,SigmaV,lambda] = pca(oriCor);
    krot=2; % Get axies
    [Urot,T] = rotatefactors(U(:,1:krot)); % get the rotation factor

    % Rotate the data
    Vrot=oriCor*Urot;  
    minVrot = min(Vrot); % Get minmun of the dataset
    sVrot = Vrot-minVrot; % shift all points to make the points attach to the x and y axies
end
%% end

out.N1=sum(inmask1);
out.N2=sum(inmask2);


%% align sites (currently it is quite rough)
idxPro1 = length(xm1);
allLocs = [xm1 ym1; xm2 ym2];
ori = [median(allLocs(:,1),1) prctile(allLocs(:,2),95,1)];
newOri = allLocs-ori;
out.z1 = newOri(1:idxPro1,2); out.z2 = newOri((idxPro1+1):end,2);
out.x1 = newOri(1:idxPro1,1); out.x2 = newOri((idxPro1+1):end,1);
out.zmed1=median(ym1);out.zmed2=median(ym2);
out.zmedd=out.zmed1-out.zmed2;
out.xmed1=median(xm1);out.xmed2=median(xm2);

kxm1 = fitdist(xm1, 'Kernel');
kxm1 = fitdist(xm1, 'Kernel', 'Width', kxm1.BandWidth * 0.7);
pdfKxm1 = pdf(kxm1, floor(min(xm1)):ceil(max(xm1)));
%figure(5556)
%plot(pdfKxm1)

kxm2 = fitdist(xm2, 'Kernel');
kxm2 = fitdist(xm2, 'Kernel', 'Width', kxm2.BandWidth * 0.7);
pdfKxm2 = pdf(kxm2, floor(min(xm2)):ceil(max(xm2)));
lmKxm2 = islocalmax(pdfKxm2);
dia = diff(find(lmKxm2));
out.dPeak = dia;
%figure(5555)
%plot(pdfKxm2)

kym1 = fitdist(ym1, 'Kernel');
pdfKym1 = pdf(kym1, floor(min(ym1)):ceil(max(ym1)));



out.zm1=mean(ym1);out.zm2=mean(ym2);
out.zmd=out.zmed1-out.zmed2;
out.xm1=mean(xm1);out.xm2=mean(xm2);
sr=p.se_siteroi;
n=-sr/2+p.se_sitepixelsize/2:p.se_sitepixelsize:sr/2;
[X,Y]=meshgrid(n,n);

out.maskmx=sum(sum(mask.*X))/sum(mask(:));
out.maskmz=sum(sum(mask.*Y))/sum(mask(:));

%quantiles
pq=[0.05 0.25 0.5 0.75 0.95];
if ~isempty(xm1)
    out.qx1=myquantile(xm1,pq);
    out.qz1=myquantile(ym1,pq);
    out.x1q13=out.qx1(4)-out.qx1(2);
    out.dz1=out.qz1(5)-out.qz1(1);
    out.z1q13=out.qz1(4)-out.qz1(2);
else
    out.qx1=[0 0 0];
    out.qz1=[0 0 0];
    out.dz1=0;
end
if ~isempty(xm2)
    out.qx2=myquantile(xm2,pq);
    out.x2q13=out.qx2(4)-out.qx2(2);
    out.qz2=myquantile(ym2,pq);
    out.dz2=out.qz2(5)-out.qz2(1);
    out.z2q13=out.qz2(4)-out.qz2(2);
else
    out.qx2=[0 0 0];
    out.qz2=[0 0 0];
    out.dz2=0;
end
    out.centdist=out.qz2(2)-out.qz1(2);


%plot all
img1=makeimage(p,locs1.xnmrot,locs1.ynmrot,locs1.locprecnm*p.layer1_.gaussfac);
img2=makeimage(p,locs2.xnmrot,locs2.ynmrot,locs2.locprecnm*p.layer2_.gaussfac);
p.saturation=2;
 hf=obj.setoutput('filtered');
imf=mask2im(p,imf1,imf2,mask1,mask2,mask12,mask);
imagesc([-roisize(1) roisize(1)],[-roisize(2) roisize(2)],imf,'Parent',hf);

 hb=obj.setoutput('image');
imc=mask2im(p,img1,img2,mask1,mask2,mask12,mask);
imagesc([-roisize(1) roisize(1)],[-roisize(2) roisize(2)],imc,'Parent',hb);
%  h2=obj.setoutput('xx');
 h2=hb;
 h2.NextPlot='add';
plot(locs1.xnmrot(inmask1),locs1.ynmrot(inmask1),'bo','Parent',h2);
plot(locs1.xnmrot(~inmask1),locs1.ynmrot(~inmask1),'b.','Parent',h2);
plot(locs2.xnmrot(inmask2),locs2.ynmrot(inmask2),'mo','Parent',h2);
plot(locs2.xnmrot(~inmask2),locs2.ynmrot(~inmask2),'m.','Parent',h2);

plot(out.xmed1,out.zmed1,'ko',out.xmed1,out.zmed1,'wx','Parent',h2);
plot(out.xmed2,out.zmed2,'ko',out.xmed2,out.zmed2,'wx','Parent',h2);
plot(out.xm1,out.zm1,'kd',out.xm1,out.zm1,'w*','Parent',h2);
plot(out.xm2,out.zm2,'kd',out.xm2,out.zm2,'w*','Parent',h2);

plot(out.maskmx,out.maskmz,'wo','Parent',h2);
h2.NextPlot='replace';
axis(h2,'ij');

%make shifted images from filtered localizations
% if ~isnan(out.xmed2)
%     xref=out.xmed2;
% else
%     xref=0;
% end
% if isfield(out,'qz2')
%     yref=out.qz2(3);
% else yref=0;
% end
% 
% x1=locs1.xnmrot(inmask1)-xref;y1=locs1.ynmrot(inmask1)-yref;
% x2=locs2.xnmrot(inmask2)-xref;y2=locs2.ynmrot(inmask2)-yref;
% img1=makeimage(p,x1,y1,locs1.locprecnm(inmask1)*p.layer1_.gaussfac,[],2*roisize);
% img2=makeimage(p,x2,y2,locs2.locprecnm(inmask2)*p.layer2_.gaussfac,[],2*roisize);
% s=size(img1);
% imgt=zeros(s(1),s(2),3);
% imgtp=imgt;
% imgt(:,:,1)=img1;
% imgt(:,:,2)=img2;
% 
% imgtp(:,:,1)=img1/max(img1(:));
% imgtp(:,:,2)=img2/max(img2(:));
% 
% 
%  hc=obj.setoutput('image_shifted');
%  imagesc([-roisize(1) roisize(1)],[-roisize(2) roisize(2)],imgtp,'Parent',hc);
%  out.imageshift=imgt;
end
function imc=mask2im(p,img1,img2,mask1,mask2,mask12,mask)
sim=size(mask1);
imc=zeros(sim(1),sim(2),3);
img1n=img1/max(img1(:))*p.saturation;
img2n=img2/max(img2(:))*p.saturation;
img1n(img1n>1)=1;img2n(img2n>1)=1;

p1=bwperim(mask1);p2=bwperim(mask2);pm=bwperim(mask);p12=bwperim(mask12);
img1n(p1)=1;img2n(p2)=1;
img1n(pm)=1;img2n(pm)=1;
imc(:,:,1)=img1n;
imc(:,:,2)=img2n;
imc(:,:,3)=pm|p12;
imc=imc/max(imc(:));
end
function ind=inmask(p,locs,mask)
roisize=ones(2,1)*p.se_siteroi/2;
pixels=p.se_sitepixelsize;
ind=withinmask(mask',(locs.xnmrot+roisize(1))/pixels,(locs.ynmrot+roisize(2))/pixels);
end
function [mask,mask1,mask2,mask12,img1,img2, PCArot1, PCArot2]=makemaskall(p,obj)
roisize=ones(2,1)*p.se_siteroi/2;
roisizeh=p.se_siteroi/2;
pixels=p.se_sitepixelsize;

locs1=obj.getLocs({'xnm','ynm','locprecnm','xnmrot','ynmrot'},'layer',1,'size',roisizeh,'grouping','ungrouped');
locs2=obj.getLocs({'xnm','ynm','locprecnm','xnmrot','ynmrot'},'layer',2,'size',roisizeh,'grouping','ungrouped');

% Yu-le added this part
if 0
    numLocs1 = length(locs1.xnmrot);
    numLocs2 = length(locs2.xnmrot);
    
    oriCor = [locs1.xnmrot locs1.ynmrot;locs2.xnmrot locs2.ynmrot];

    % Do PCA together
    [U,SigmaV,lambda] = pca(oriCor);

    krot=2; % Get axies
    [Urot,T] = rotatefactors(U(:,1:krot)); % get the rotation factor

    % Rotate the data
    Vrot=oriCor*Urot;  
    minVrot = min(Vrot); % Get minmun of the dataset
    sVrot = Vrot-minVrot; % shift all points to make the points attach to the x and y axies

    theta = -pi/2;
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];

    PCArot1 = sVrot(1:numLocs1,:);
    PCArot2 = sVrot((numLocs1+1):(numLocs1+numLocs2),:);

    PCArot1 = PCArot1 * R;
    PCArot2 = PCArot2 * R;

    xm1=PCArot1(:,1);ym1=PCArot1(:,2); % Get coordinates in ROI
    xm1m = median(xm1); ym1m = median(ym1);
    xm2=PCArot2(:,1);ym2=PCArot2(:,2); % The same
    xm2m = median(xm2); ym2m = median(ym2);
end

xm1=locs1.xnmrot;ym1=locs1.ynmrot; % Get coordinates in ROI
xm2=locs2.xnmrot;ym2=locs2.ynmrot; % The same

img1=makeimage(p,xm1,ym1,p.sigma1(1),p.sigma1(2)); % convert coordinates into pixel matrix
img2=makeimage(p,xm2,ym2,p.sigma2(1),p.sigma2(2)); % The same
maxone1=1/(2*pi*p.sigma1(1)*p.sigma1(2))*pixels^2;
maxone2=1/(2*pi*p.sigma2(1)*p.sigma2(2))*pixels^2;

cutoff=p.cutofffactor;
if length(cutoff)==1;
    cutoff(2)=cutoff(1);
end
mask1=makemask(p,img1,maxone1*cutoff(1));



mask=false(size(mask1));
if p.segmentl1
    mask=mask|mask1;
end
if p.segmentl2
    mask2=makemask(p,img2,maxone2*cutoff(2));
    mask=mask|mask2;
else
    mask2=false(size(mask1));
end
if p.segmentl1l2
    mask12=makemask(p,img1/cutoff(1)/maxone1+img2/cutoff(2)/maxone2,1);
    mask=mask|mask12;
else
    mask12=false(size(mask1));
end
end

function [mask,im1,cutoff]=makemask(p,im1,maxone)
% 
% im1=makeimage(p,xm1,ym1,s(1),s(2));

p.take2factor=1.5;
cutoff= maxone;% *p.cutofffactor;
im1bw=im1>cutoff;

% if two largest segments are similar in size
im1bwa1=bwareafilt(im1bw,1); % connected objects 1
im1bwa2=bwareafilt(im1bw,2); % connected objects 1+2
if sum(im1bwa2(:))>p.take2factor*sum(im1bwa1(:))
    im1bwa=im1bwa2;
else
    im1bwa=im1bwa1;
end

sel=strel('disk',p.dilation); % create a basic structure

im1bwa=imdilate(im1bwa,sel);
mask=imfill(im1bwa,'holes');
end

function im=makeimage(p,xm,ym,sx,sy,roisize)
if nargin<6
    roisize=p.se_siteroi/2;
end
if nargin<5||isempty(sy)
    sy=sx;
end
if length(sx)==1||length(sy)==1
    sx=sx+0*xm; sy=sy+0*xm;
end

pixels=p.se_sitepixelsize;
 range=[-roisize(1) roisize(1)];
 posf.x=xm;posf.y=ym;posf.sx=sx;posf.sy=sy;
im=double(gaussrender_ellipt(posf,range, range, pixels, pixels));
end
function [z, bw] = getKernelMatrix(matrixSize, sVrot, varargin) 
    [meshx, meshy] = meshgrid(0:matrixSize(1), 0:matrixSize(2)); % make full grid
    q = [meshx(:), meshy(:)]; % convert the grid to positions
    if length(varargin)==0
        [sVrotx,xy,bw]=ksdensity(sVrot, q); % get kernel density estimation
    else
        [sVrotx,xy,bw]=ksdensity(sVrot, q, varargin{1,1}, varargin{1,2}); % get kernel density estimation    
    end
    
    % converted into an image
    idx = sub2ind(matrixSize(end:-1:1)+1, xy(:,2)+1, xy(:,1)+1); 
    z = zeros(matrixSize(end:-1:1)+1);
    z(idx) = sVrotx;
end
function [filteredCor, kf] = kernelFeatures(parentalLocs)
    kf = [];
    [~, bw1] = getKernelMatrix(kernelMatrixSize(parentalLocs), parentalLocs(sub,:));

    z1Ori = getKernelMatrix(kernelMatrixSize(parentalLocs), parentalLocs(sub,:), 'Bandwidth', bw1*0.7);
     % offset
    Y1 = quantile(z1Ori(:),0.8);
    z1 = z1Ori;
    z1(z1Ori<=Y1)=0;
    filteredCor = corFiltedByMask(z1, parentalLocs(sub,:));
    [~, bw1] = getKernelMatrix(kernelMatrixSize(parentalLocs), filteredCor);
    z1 = getKernelMatrix(kernelMatrixSize(parentalLocs), filteredCor, 'Bandwidth', bw1*0.7);
    
    BW = imregionalmax(z1);
    [i,j]=find(BW>=1);

    peakVal = z1(sub2ind(size(BW), i, j));
    lessThanf = find(max(peakVal)./peakVal<=1.2);
    peakLoc = [j i];
    peakLoc = peakLoc(lessThanf,:);
    numPeak = length(lessThanf);


    Q1 = quantile(filteredCor,0.25);
    Q2 = quantile(filteredCor,0.5);
    Q3 = quantile(filteredCor,0.75);
    DQ = Q3-Q1;
    shiftFactor = (Q2-Q1)./(Q3-Q2);
    
    pd = pdist(peakLoc, 'euclidean');
    
    kf.Q1 = Q1;
    kf.Q2 = Q2;
    kf.Q3 = Q3;
    kf.DQ = DQ;
    kf.shiftFactor = shiftFactor;
    kf.numPeak = numPeak;
    kf.peakLoc = peakLoc;
    kf.peakDist = pd;
end