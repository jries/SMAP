function [P,CRLB,LogL,color,LLsecond]=mleFit_LM_global(varargin)
% varargin:
%imstack, fittype, sharedflag, iterations, spline coefficients, channelshift,
%cmos varmap, silent,zstart,photratios

% 1. imagestack (single)
% 2. fittype:
%   1: Gauss
%   2: spline

% 3. shared parameters: simple vector (same for all localizations).

% optional:
% 4. iterations (default=50)
% 5. paramters for fitters:
%   cspline: cspline coefficients (single)
%   cspline for 2D PSF: as 5, but two fits to break asymmetry. cspline coefficients (single)
% 6. shifts dT (2xchannel matrix, only dx and dy for each channel)
% 7. silent (suppress output if 1)
% 8. scmoS varmap
% 9. z start parameter (more than one: return solution with maximum
%       LIkelihood). Units: distance of stack calibration, center based.
%       This is a simple vector, same for all localizations.
% 10. Intenstiy ratios to test for multi-color SMLM, don't link if empty

%Output:
%P (fittype): multiple copies for multiple channels
%1. X, Y, Photons, Background, PSFxy, Iterations
%2. X, Y, Photons, Background, Z, Iterations

%CRLB: cramer-rao lower bounds, as in P
%LogL: log-likelihood.

%Only for fitmode 6: P1 etc: results with z-startparameter<0, P2 etc:
%results with z-startparameter>0

color=[];


imstack=single(varargin{1});
Nfits = size(imstack,3);
fittype = int32(varargin{2});
sharedin=int32(varargin{3});
shared = repmat(int32(sharedin),[1 Nfits]); %for GPU fitter: avoid access to same variables from multiple threads

iterations=int32(varargin{4});
splinecoeff=single(varargin{5});
channelshift=single(varargin{6});
coeffsize=size(splinecoeff);

if nargin<7 || isempty(varargin{8}) %silent
    silent=0;
else
    silent=single(varargin{8});
end

if nargin<8 || isempty(varargin{7}) %varmap
    varmap=single(varargin{7});
else
    varmap=0;
end

if nargin<9||isempty(varargin{9}) %zstart
    zstart=single(coeffsize(3)/2);
else
    zstart=single(coeffsize(3)/2+varargin{9});
end

if nargin<10 || isempty(varargin{10}) %fixed photon ratios
    PhotonRatios=1;
    photonratiofixed=false;
else
    PhotonRatios=varargin{10};
    shared(4,:)=1; %if we fit with ratios, N needs to be linked
    iterationsin=iterations;
    iterations=15; %used for testing different colors
    photonratiofixed=true;
end
dT=ratiochannelshift(channelshift,1); %for test and used if no photon ratio given

%determine of it runs on GPU, otherwise use CPU as default
persistent fitter
allfitters={@GPUmleFit_LM_MultiChannel,@CPUmleFit_LM_MultiChannel};
if isempty(fitter)
    for k=1:length(allfitters)
        try
            allfitters{k}(imstack,fittype,shared,iterations,splinecoeff,dT);
            fitter=k;
            break
        catch err
            % fitter did not work
        end
    end
    disp(['using: ' char(allfitters{fitter})]);
end

if isempty(fitter)
    error('No working GPU or CPU fitter found.')
end
fitterh=allfitters{fitter};

%do fitting
first=true;
LogLp=ones(size(imstack,3),length(PhotonRatios));
for phot=1:length(PhotonRatios)
    if photonratiofixed
        dT=ratiochannelshift(channelshift,PhotonRatios(phot));
    end
    for k=1:length(zstart)
         zstarth=repmat(zstart(k),1,Nfits);
         [Ph,CRLBh,LogLh]=fitterh(imstack,fittype,shared,(iterations),splinecoeff,dT,varmap,silent,zstarth);
        if first %first round
            P=Ph;
            CRLB=CRLBh;
            LogL=LogLh;
            first=false;
            color=zeros(size(LogLh))+phot;
            zpos=ones(size(LogLh));
        else
            indbetter=LogLh-LogL>0; %copy only everything if LogLh increases by more than rounding error.
            P(indbetter,:)=Ph(indbetter,:);
            CRLB(indbetter,:)=CRLBh(indbetter,:);
            LogL(indbetter)=LogLh(indbetter);
            color(indbetter)=phot;
            zpos(indbetter)=zstart(k);
        end
    end
    LogLp(:,phot)=LogLh;
end

if photonratiofixed
    zstarth=P(:,3);
    dT(4,4,:)=PhotonRatios(color);
    [P,CRLB,LogL]=fitterh(imstack,fittype,shared,iterationsin,splinecoeff,dT,varmap,silent,zstarth);
    ml=max(LogLp,[],2);
    LLsecond=quantile(ml./LogLp,1-1/length(PhotonRatios),2); 
else
    color=[];LLsecond=[];
end
end

function  dT=ratiochannelshift(channelshift,PhotonRatio)
dS1 = [1, 1 ;1, 1 ;1, 1;1, PhotonRatio;1, 1];
npar=5;
Nfits=size(channelshift,3);
noChannels=size(channelshift,2);
temp1 = repmat(dS1,[1 1 Nfits]);
dT = zeros(npar,noChannels*2,Nfits,'single');
dT(1:2,1:noChannels,:)=channelshift(1:2,:,:);
dT(:,noChannels+1:2*noChannels,:)=temp1;
end
 
