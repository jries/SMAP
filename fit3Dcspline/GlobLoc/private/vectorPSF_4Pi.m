function [allPSFs,allPSFs_upper, allPSFs_lower,Waberration_upper,Waberration_lower] = vectorPSF_4Pi(parameters)
% parameters: NA, refractive indices of medium, cover slip, immersion fluid,
% wavelength (in nm), sampling in pupil
%optics para
% add norm function comparecd to v3 #need test
% add pixelX Y

NA = parameters.NA;
refmed = parameters.refmed;
refcov = parameters.refcov;
refimm = parameters.refimm;
lambda = parameters.lambda;

%image para
xemit_upper = parameters.xemit;
yemit_upper = parameters.yemit;
zemit_upper = parameters.zemit;
objStage_upper=parameters.objStage;

objStage0_upper=parameters.objStage0_upper;
zemit0_upper = parameters.zemit0_upper;

xemit_lower = parameters.xemit+parameters.offset(1);
yemit_lower = parameters.yemit+parameters.offset(2);
zemit_lower = -1.*parameters.zemit;
objStage_lower=-1*parameters.objStage;

objStage0_lower=parameters.objStage0_lower;
zemit0_lower = parameters.zemit0_lower;

phaseshift = parameters.phaseshift;

% objStage0 = parameters.objStage0;
Npupil = parameters.Npupil;
sizeX = parameters.sizeX;
sizeY = parameters.sizeY;
sizeZ = parameters.Nmol;
pixelSizeX = parameters.pixelSizeX;
pixelSizeY = parameters.pixelSizeY;
xrange = pixelSizeX*sizeX/2;
yrange = pixelSizeY*sizeY/2;

% pupil radius (in diffraction units) and pupil coordinate sampling
PupilSize = 1.0;
DxyPupil = 2*PupilSize/Npupil;
XYPupil = -PupilSize+DxyPupil/2:DxyPupil:PupilSize;
[YPupil,XPupil] = meshgrid(XYPupil,XYPupil);

% calculation of relevant Fresnel-coefficients for the interfaces
% between the medium and the cover slip and between the cover slip
% and the immersion fluid
CosThetaMed = sqrt(1-(XPupil.^2+YPupil.^2)*NA^2/refmed^2);
CosThetaCov = sqrt(1-(XPupil.^2+YPupil.^2)*NA^2/refcov^2);
CosThetaImm = sqrt(1-(XPupil.^2+YPupil.^2)*NA^2/refimm^2);
%need check again, compared to original equation, FresnelPmedcov is
%multipiled by refmed
FresnelPmedcov = 2*refmed*CosThetaMed./(refmed*CosThetaCov+refcov*CosThetaMed);
FresnelSmedcov = 2*refmed*CosThetaMed./(refmed*CosThetaMed+refcov*CosThetaCov);
FresnelPcovimm = 2*refcov*CosThetaCov./(refcov*CosThetaImm+refimm*CosThetaCov);
FresnelScovimm = 2*refcov*CosThetaCov./(refcov*CosThetaCov+refimm*CosThetaImm);
FresnelP = FresnelPmedcov.*FresnelPcovimm;
FresnelS = FresnelSmedcov.*FresnelScovimm;

% Apoidization for sine condition
apoid = sqrt(CosThetaImm)./CosThetaMed;
% definition aperture
ApertureMask = double((XPupil.^2+YPupil.^2)<1.0);
Amplitude = ApertureMask.*apoid;

% setting of vectorial functions
Phi = atan2(YPupil,XPupil);
CosPhi = cos(Phi);
SinPhi = sin(Phi);
CosTheta = sqrt(1-(XPupil.^2+YPupil.^2)*NA^2/refmed^2);
SinTheta = sqrt(1-CosTheta.^2);

pvec{1} = FresnelP.*CosTheta.*CosPhi;
pvec{2} = FresnelP.*CosTheta.*SinPhi;
pvec{3} = -FresnelP.*SinTheta;
svec{1} = -FresnelS.*SinPhi;
svec{2} = FresnelS.*CosPhi;
svec{3} = 0;

PolarizationVector = cell(2,3);
for jtel = 1:3
    PolarizationVector{1,jtel} = CosPhi.*pvec{jtel}-SinPhi.*svec{jtel};
    PolarizationVector{2,jtel} = SinPhi.*pvec{jtel}+CosPhi.*svec{jtel};
end

wavevector = cell(1,3);
wavevector{1} = (2*pi*NA/lambda)*XPupil;
wavevector{2} = (2*pi*NA/lambda)*YPupil;
wavevector{3} = (2*pi*refimm/lambda)*CosThetaImm;
wavevectorzmed = (2*pi*refmed/lambda)*CosThetaMed;

% calculation aberration function
Waberration_lower = zeros(size(XPupil));
Waberration_upper = zeros(size(XPupil));
orders_lower = parameters.aberrations(:,1:2,1);
orders_upper = parameters.aberrations(:,1:2,2);
zernikecoefs_lower = squeeze(parameters.aberrations(:,3,1));
zernikecoefs_upper = squeeze(parameters.aberrations(:,3,2));

normfac_lower = sqrt(2*(orders_lower(:,1)+1)./(1+double(orders_lower(:,2)==0)));
normfac_upper = sqrt(2*(orders_upper(:,1)+1)./(1+double(orders_upper(:,2)==0)));

zernikecoefs_lower = normfac_lower.*zernikecoefs_lower;
zernikecoefs_upper = normfac_upper.*zernikecoefs_upper;

allzernikes_lower = get_zernikefunctions(orders_lower,XPupil,YPupil);
for j = 1:numel(zernikecoefs_lower)
    Waberration_lower = Waberration_lower+zernikecoefs_lower(j)*squeeze(allzernikes_lower(j,:,:));
end


allzernikes_upper = get_zernikefunctions(orders_upper,XPupil,YPupil);
for j = 1:numel(zernikecoefs_upper)
    Waberration_upper = Waberration_upper+zernikecoefs_upper(j)*squeeze(allzernikes_upper(j,:,:));
end

% pupil and image size (in diffraction units)
% PupilSize = NA/lambda;
ImageSizex = xrange*NA/lambda;
ImageSizey = yrange*NA/lambda;

% calculate auxiliary vectors for chirpz
[Ax,Bx,Dx] = prechirpz(PupilSize,ImageSizex,Npupil,sizeX);
[Ay,By,Dy] = prechirpz(PupilSize,ImageSizey,Npupil,sizeY);
allPSFs = zeros(sizeX,sizeY,sizeZ,length(phaseshift));
allPSFs_upper = zeros(sizeX,sizeY,sizeZ,length(phaseshift));
allPSFs_lower = zeros(sizeX,sizeY,sizeZ,length(phaseshift));
for ii = 1:length(phaseshift)
    pphase = phaseshift(ii);
    PhaseFactor_lower = exp(2*pi*1i*Waberration_lower/lambda).*exp(-1i*pphase); % compare to the flip version, we have to use neg sign here
    PhaseFactor_upper = exp(2*pi*1i*Waberration_upper/lambda);
    
     % compute pupil matrix
    PupilMatrix_lower = cell(2,3);
    for itel = 1:2
        for jtel = 1:3
            PupilMatrix_lower{itel,jtel} = Amplitude.*PhaseFactor_lower.*PolarizationVector{itel,jtel};
        end
    end
    
    PupilMatrix_upper = cell(2,3);
    for itel = 1:2
        for jtel = 1:3
            PupilMatrix_upper{itel,jtel} = Amplitude.*PhaseFactor_upper.*PolarizationVector{itel,jtel};
        end
    end
    
    
    
    FieldMatrix_lower = cell(2,3,sizeZ);
    FieldMatrix_upper = cell(2,3,sizeZ);
    

    
    for jz = 1:sizeZ
        % xyz induced phase contribution
        if zemit_lower(jz)+zemit0_lower>=0
            Wxyz_lower= (-1*xemit_lower(jz))*wavevector{1}+(-1*yemit_lower(jz))*wavevector{2}+(zemit_lower(jz)+zemit0_lower)*wavevectorzmed;
            %   zemitrun = objStage(jz);
            PositionPhaseMask_lower = exp(1i*(Wxyz_lower+(objStage_lower(jz)+objStage0_lower)*wavevector{3}));
        else
            Wxyz_lower= (-1*xemit_lower(jz))*wavevector{1}+(-1*yemit_lower(jz))*wavevector{2};
            %   zemitrun = objStage(jz);
            PositionPhaseMask_lower = exp(1i*(Wxyz_lower+(objStage_lower(jz)+objStage0_lower+zemit_lower(jz)+zemit0_lower)*wavevector{3}));
        end
        
        for itel = 1:2
            for jtel = 1:3
                
                % pupil functions and FT to matrix elements
                PupilFunction_lower =PositionPhaseMask_lower.*PupilMatrix_lower{itel,jtel};
                IntermediateImage_lower = transpose(cztfunc(PupilFunction_lower,Ay,By,Dy));
                FieldMatrix_lower{itel,jtel,jz} = transpose(cztfunc(IntermediateImage_lower,Ax,Bx,Dx));
            end
        end
    end
    
    for jz = 1:sizeZ
        % xyz induced phase contribution
        if zemit_upper(jz)+zemit0_upper>=0
            Wxyz_upper= (-1*xemit_upper(jz))*wavevector{1}+(-1*yemit_upper(jz))*wavevector{2}+(zemit_upper(jz)+zemit0_upper)*wavevectorzmed;
            %   zemitrun = objStage(jz);
            PositionPhaseMask_upper = exp(1i*(Wxyz_upper+(objStage_upper(jz)+objStage0_upper)*wavevector{3}));
        else
            Wxyz_upper= (-1*xemit_upper(jz))*wavevector{1}+(-1*yemit_upper(jz))*wavevector{2};
            %   zemitrun = objStage(jz);
            PositionPhaseMask_upper = exp(1i*(Wxyz_upper+(objStage_upper(jz)+objStage0_upper+zemit_upper(jz)+zemit0_upper)*wavevector{3}));
        end
        
        for itel = 1:2
            for jtel = 1:3
                
                % pupil functions and FT to matrix elements
                PupilFunction_upper =PositionPhaseMask_upper.*PupilMatrix_upper{itel,jtel};
                IntermediateImage_upper = transpose(cztfunc(PupilFunction_upper,Ay,By,Dy));
                FieldMatrix_upper{itel,jtel,jz} = transpose(cztfunc(IntermediateImage_upper,Ax,Bx,Dx));
            end
        end
    end
    
    
    FieldMatrixAll = FieldMatrix_upper;
    for jz= 1:size(FieldMatrixAll,3)
        for i = 1:2
            for j= 1:3
                FieldMatrixAll{i,j,jz} = FieldMatrix_upper{i,j,jz}+FieldMatrix_lower{i,j,jz};
            end
        end
    end
    
    
    
    %calculates the free dipole PSF given the field matrix.
    PSFs_upper = zeros(sizeX,sizeY,sizeZ);
    for jz = 1:sizeZ
        for jtel = 1:3
            for itel = 1:2
                PSFs_upper(:,:,jz) = PSFs_upper(:,:,jz) + (1/3)*abs(FieldMatrix_upper{itel,jtel,jz}).^2;
            end
        end
    end
    
    PSFs_lower = zeros(sizeX,sizeY,sizeZ);
    for jz = 1:sizeZ
        for jtel = 1:3
            for itel = 1:2
                PSFs_lower(:,:,jz) = PSFs_lower(:,:,jz) + (1/3)*abs(FieldMatrix_lower{itel,jtel,jz}).^2;
            end
        end
    end
    
    PSFs = zeros(sizeX,sizeY,sizeZ);
    for jz = 1:sizeZ
        for jtel = 1:3
            for itel = 1:2
                PSFs(:,:,jz) = PSFs(:,:,jz) + (1/3)*abs(FieldMatrixAll{itel,jtel,jz}).^2;
            end
        end
    end
    
    
    
    
    % calculate intensity normalization function using the PSFs at focus
    % position without any aberration. It might not work when sizex and sizeY
    % are very small
    FieldMatrix = cell(2,3);
    for itel = 1:2
        for jtel = 1:3
            PupilFunction = Amplitude.*PolarizationVector{itel,jtel};
            IntermediateImage = transpose(cztfunc(PupilFunction,Ay,By,Dy));
            FieldMatrix{itel,jtel} = transpose(cztfunc(IntermediateImage,Ax,Bx,Dx));
        end
    end
    
    intFocus = zeros(sizeX,sizeY);
    for jtel = 1:3
        for itel = 1:2
            intFocus = intFocus + (1/3)*abs(FieldMatrix{itel,jtel}).^2;
        end
    end
    
    normIntensity = sum(intFocus(:));
    
    PSFs = PSFs./normIntensity./2;
    
    PSFs_upper = PSFs_upper./normIntensity;
    
    PSFs_lower = PSFs_lower./normIntensity;
    
    
%     h = fspecial('gaussian',5,0.5);
%     PSFs = convn(PSFs,h,'same');
    
    allPSFs(:,:,:,ii) = PSFs;
    allPSFs_upper(:,:,:,ii) = PSFs_upper;
    allPSFs_lower(:,:,:,ii) = PSFs_lower;
    
    
end






