classdef CME2CSide_yule<interfaces.SEEvaluationProcessor
    properties
        boundary
    end
    methods
        function obj=CME2CSide_yule(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj, inp)
            roisize=obj.getPar('se_siteroi');
            locsL1=obj.getLocs({'locprecnm','PSFxnm','xnm','ynm','frame','channel', 'gaussfac'},'layer',1);
            locsL2=obj.getLocs({'locprecnm','PSFxnm','xnm','ynm','frame','channel', 'gaussfac'},'layer',2);
            oriCor = [locsL1.xnm locsL1.ynm; locsL2.xnm locsL2.ynm];
                       
            
            
            % Do PCA together
            [U,SigmaV,lambda] = pca(oriCor);
            
            
            krot=2; % Get axies
            [Urot,T] = rotatefactors(U(:,1:krot)); % get the rotation factor
            
            % Rotate the data
            Vrot=oriCor*Urot;  
            
            minVrot = min(Vrot); % Get minmun of the dataset

            sVrot = Vrot-minVrot; % shift all points to make the points attach to the x and y axies
            
             % Get the kernel matrix
            img = makeimage(inp, sVrot(:,1), sVrot(:,2), 15, 5, 500);
            
            pixels=inp.se_sitepixelsize;
            maxone1=1/(2*pi*15*5)*pixels^2;

            cutoff=1;
            if length(cutoff)==1
                cutoff(2)=cutoff(1);
            end
            inp.dilation = 2;
            mask1=makemask(inp,img,maxone1*cutoff(1));            
            
            
            %[inLoc1, kf1] = kernelFeatures(sVrot, locs.channel==1);
            %[inLoc2, kf2] = kernelFeatures(sVrot, locs.channel==2);
            

            
            
            figure(556)
            clf
            imagesc(z1)
            hold on
            
            
            figure(559)
            clf
            [C1,H1] = contour(z1, 200);
            hold on
            L1 = 70;
            contour(z1,[H1.LevelList(L1), H1.LevelList(L1)],'LineWidth',2 );
            
          
            
            scatter(j,i)
            scatter(sVrot(:,1), sVrot(:,2), 'MarkerFaceColor',[0 0 0])
            
            
            
            
             z2 = getKernelMatrix(sVrot(locs.channel==2,:), 'Bandwidth', bw2*0.7);
            
            Y2 = quantile(z2(:),0.95);
            z2(z2<=Y2)=0;
            
            
            
            

            
            figure(555)
            clf
            imagesc(z2)
            hold on
            BW = imregionalmax(z2);
            imagesc(BW)
            contour(z2, 200)
            contour(z1, 200)
            
            scatter(sVrot(:,1), sVrot(:,2), 'MarkerFaceColor',[0 0 0])
            
            
            
            
            % offset
            Y = mean(z(:))*1.5;
            z(z<=Y)=0;

            figure(556)
            clf
            BW = imregionalmax(z);
            hold on
            imagesc(z)
            scatter(sVrot(:,1), sVrot(:,2), 'MarkerFaceColor',[0 0 0])

        end

        function pard=guidef(obj)
            pard=guidef;
        end
    end

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
            
function pard=guidef
pard.t_gridRange.object=struct('Style','text','String','Bin size of the grids');
pard.t_gridRange.position=[1,1];
pard.t_gridRange.Width=2;

pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi','layer1_','layer2_','se_sitepixelsize'};
end

function filtered = corFiltedByMask(mask, cor)
    sub = find(mask);
    ccor = ceil(cor);
    q = sub2ind(size(mask), ccor(:,2)'+1, ccor(:,1)'+1);
    filtered = ismember(q,sub);
    filtered = cor(filtered,:);
end

function maxSVrot = kernelMatrixSize(sVrot)
    maxSVrot = max(sVrot);
    maxSVrot = ceil(maxSVrot); % get maximun of x and y position
end

function [filteredCor, kf] = kernelFeatures(parentalLocs, sub)
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

function [mask,im1,cutoff]=makemask(p,im1,maxone)
% 
% im1=makeimage(p,xm1,ym1,s(1),s(2));

p.take2factor=1.5;
cutoff= maxone;% *p.cutofffactor;
im1bw=im1>cutoff;

% if two largest segments are similar in size
im1bwa1=bwareafilt(im1bw,1);
im1bwa2=bwareafilt(im1bw,2);
if sum(im1bwa2(:))>p.take2factor*sum(im1bwa1(:))
    im1bwa=im1bwa2;
else
    im1bwa=im1bwa1;
end

sel=strel('disk',p.dilation);

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
function ind=inmask(p,locs,mask)
roisize=ones(2,1)*p.se_siteroi/2;
pixels=p.se_sitepixelsize;
ind=withinmask(mask',(locs.xnmrot+roisize(1))/pixels,(locs.ynmrot+roisize(2))/pixels);
end