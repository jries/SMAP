%%
% Output
% Results = [frame,,xPos,yPos, amplitude, sigma,background,resnorm/sum(ccd_image(:))]
%%
function [Results] = LS_fittingVer2(filename,num_frame,dummy_frame,ADU,baseline,thresh,debug)
Results = [];
num_frame = min(num_frame,length(imfinfo(filename)));
fprintf('running \n');
parfor ii = 1:num_frame
    fprintf('running... %d frame \n',ii);
    low_img = imread(filename,ii+dummy_frame);
    low_img = double(ADU*(low_img-baseline));
    if debug > 0
        figure(1); imagesc(low_img); colormap(hot); axis image; colorbar;
    end 
    [LS_result]= Gaussfit(low_img,thresh);
    if~isempty(LS_result)
        Result_temp = zeros(size(LS_result,1),7);
        Result_temp(:,1) = ii;
        Result_temp(:,2:3) = LS_result(:,3:4)-0.5;
        Result_temp(:,4) = LS_result(:,1);
        Result_temp(:,5) = LS_result(:,2);
        Result_temp(:,6) = LS_result(:,5);
        Result_temp(:,7) = LS_result(:,6);
        Results = [Results;Result_temp];
    end
end

% fitting the data by sum of two gaussians.
function [fitresults]= Gaussfit(image,threshold)

epsilon1 =  10^-3; 
epsilon2 =  10^-3;
maxIter  = 50;

curvefitoptions = optimset( 'lsqcurvefit');
curvefitoptions = optimset( curvefitoptions,'Jacobian' ,'off','Display', 'off',  'TolX', epsilon2, 'TolFun', epsilon1,'MaxPCGIter',1,'MaxIter',maxIter);
fitresults=[];
temp= calcstorm(image,threshold,curvefitoptions); % do the storm fit %MATLAB VERSION
fitresults=[fitresults;temp];


function fitresults = calcstorm(image,threshold,curvefitoptions) % do the storm calculation

windowSize= 3;
[sizey sizex] = size(image);

points = findpeaksMod(image, threshold);
if isempty(points)==0
    nrspots=length(points(:,1));
    
    fitresults=zeros(nrspots,6);
    
    for i=1:nrspots,
        X0=points(i,1);
        Y0=points(i,2);
        
        %round X0, Y0 to use as matrix locations
        X0_int = round(X0);
        Y0_int = round(Y0);
        
        xstart =  X0_int-windowSize;
        xfinish = X0_int+windowSize;
        ystart =  Y0_int-windowSize;
        yfinish = Y0_int+windowSize;
        
        if (xstart>1) &&  (xfinish < sizex) && (ystart>1) && (yfinish < sizey)
            
            X_POSim = X0-xstart+1;
            Y_POSim = Y0-ystart+1;
            
            img=double(image(ystart:yfinish,xstart:xfinish));
            background= min(img(:));
            brightness = sum(img(:)-background);
            widthStart = 1.5;

            initguess=double([brightness,widthStart,X_POSim,Y_POSim, background]);
            xLim = [0 sizex];
            yLim = [0 sizey];
            sigmaLim = [0.5  3];
            
            [fitParams,res] = doubleGauss2d_Fit(img,initguess ,xLim,yLim, sigmaLim,curvefitoptions);%MATLAB ONLY
                        
            % assign the data
            xPos = fitParams(3)+xstart-1;
            yPos = fitParams(4)+ystart-1;
            sigma = fitParams(2);
            background = fitParams(5);
            amplitude = fitParams(1);

            fitresults(i,:)=[amplitude,sigma,xPos,yPos,background,res/sum(img(:))];
        end
    end
    fitresults = fitresults(all(fitresults,2) >0,:);
else
    fitresults=[];
end

%----------------------------------------------------------------------------------------------
function [fitParam,res] = doubleGauss2d_Fit(inputIm,initguess ,xLim,yLim, sigmaLim,curvefitoptions)
Astart = initguess(1);
widthStart = initguess(2);
xStart = initguess(3);
yStart = initguess(4);
BGstart = initguess(5);

xMin = xLim(1);
xMax = xLim(2);
yMin = yLim(1);
yMax = yLim(2);

sigmaMin = sigmaLim(1);
sigmaMax = sigmaLim(2);
[sizey ,sizex] = size(inputIm);
[X,Y]= meshgrid(1:sizex,1:sizey);
grid = [X Y];

initGuess5Vector = [Astart widthStart  xStart yStart BGstart];

lb = [ 1        sigmaMin   xMin yMin    0  ];
ub = [65535     sigmaMax   xMax yMax  65535];

try
    [fitParam, res] = ...
        lsqcurvefit(@(x, xdata) doubleGauss2d(x, xdata), ...
        initGuess5Vector ,grid ,inputIm ,...
        lb,ub,curvefitoptions);
catch ME
    if strcmp(ME.identifier,'optim:snls:InvalidUserFunction') 
        fitParam = [0 0 0 0 0 0 0];
        res = 0;
    else
        rethrow(ME);
    end
end

if fitParam(1)< 0 
    fitParam = [0 0 0 0 0 0 0];
end


function F = doubleGauss2d(a, data)
% sum of two symmetric Gaussian functions 
% a(1) - A
% a(2) - sigma
% a(3) - Xpos
% a(4) - Ypos
% a(5) - B

[sizey sizex] = size(data);
sizex = sizex/2;

X = data(:,1:sizex);
Y = data(:,sizex+1:end);

expPart = exp(-( (X-a(3)).^2 + (Y-a(4)).^2) /(2*(a(2))^2));
F =  a(1)*expPart/(2*pi*a(2)^2) + a(5); 






%%
% below functions are related to local peak detection 

%-----------------------------------------------------------------------------------
function peaks = findpeaksMod(im, threshold)
wavelet_level = 5;
background = background_estimation(im,1,wavelet_level,'db6',5);
g = fspecial('gaussian',11,1);
im2 = imfilter(im-background,g,'symmetric','conv');
[peaks] = Peakdet2D(im2,threshold,2);

%-----------------------------------------------------------------------------------
