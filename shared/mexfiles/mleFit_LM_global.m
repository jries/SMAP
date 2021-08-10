function [P,CRLB,LogL,color,LLsecond]=mleFit_LM_global(varargin)
% varargin:
%imstack, sharedflag, iterations, spline coefficients, channelshift,
%fitmode, varmap,zstart

% 1. imagestack (single)
% 2. shared parameters

% optional:
% 3. iterations (default=50)
% 4. paramters for fitters:

%   cspline: cspline coefficients (single)
%   cspline for 2D PSF: as 5, but two fits to break asymmetry. cspline coefficients (single)
% 5. shifts dT
% not working 6. silent (suppress output if 1)
% not working7. scmoS varmap
% 6. z start parameter (more than one: return solution with maximum
%       LIkelihood). Units: distance of stack calibration, center based
% 7. Intenstiy ratios to test for multi-color SMLM, don't link if empty

%Output:
%P
%1. X, Y, Photons, Background, Iterations
%2. X, Y, Photons, Background, PSFxy, Iterations
%3. X, Y, Photons, Background, Z, Iterations
%4. X, Y, Photons, Background, PSFx, PSFy, Iterations
%5. X, Y, Photons, Background, Z, Iterations
%6. X, Y, Photons, Background, Z, Iterations
%CRLB: cramer-rao lower bounds, as in P
%LogL: log-likelihood.

%Only for fitmode 6: P1 etc: results with z-startparameter<0, P2 etc:
%results with z-startparameter>0

% [P,CRLB, LL] =CPUmleFit_LM_MultiChannel(fitstack,int32(sharedA),50,coeffh,single(dT),fitmode);
% P=[];CRLB=[];LogL=[]; %
color=[];
LLsecond=[];
%determine of it runs on GPU, otherwise use CPU as default


imstack=single(varargin{1});
shared=int32(varargin{2});
splinecoeff=single(varargin{4});
channelshift=single(varargin{5});

coeffsize=size(splinecoeff);

if nargin>2 && ~isempty(varargin{3})
    iterations=varargin{3};
else
    iterations=30;
end

if nargin<7 || isempty(varargin{7})
    PhotonRatios=1;
    dT=channelshift;
    photonratiofixed=false;
else
    PhotonRatios=varargin{7};
    shared(4,:)=1; %if we fit with ratios, N needs to be linked
    iterations=10;
    photonratiofixed=true;
    iterationsin=varargin{3};
end

if nargin<6||isempty(varargin{6})
    zstart=single(coeffsize(3)/2); %emccd
else
    zstart=single(coeffsize(3)/2+varargin{6});
end





% select fitter:    
persistent fitter
if ~photonratiofixed
    allfitters={@GPUmleFit_LM_MultiChannel_11,@CPUmleFit_LM_MultiChannel};
    allfittersnames={'GPUmleFit_LM_MultiChannel_11','CPUmleFit_LM_MultiChannel'};
    if isempty(fitter)
        for k=1:length(allfitters)
            try
    %             allfitters{k}(ones(7,7,1,2,'single'),single([1 1 1 0 0]),varargin{4},zeros(5,2));
                allfitters{k}(imstack,int32(5),shared,int32(iterations),splinecoeff,single(dT));
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
else
    fitterh=@GPUmleFit_LM_MultiChannel_vectorZ;
end

%Yiming: PhotonRatios is a vector with photon ratios. You know best how to
%pass it on to the function. If length(PhotonRatios)>1: loop over, choose
%highest LL, save color (e.g. at P(:, end+1))
%try out: test with small number of iterations, then do final fit withlarge
%iterations and zstart=zbest. For this zstart should be particle-wise

first=true;
LogLp=ones(size(imstack,3),length(PhotonRatios));
for phot=1:length(PhotonRatios)
    dT=ratiochannelshift(channelshift,PhotonRatios(phot));
%     [P,CRLB,LogL]=fitterh(imstack,int32(5),shared,int32(iterations),splinecoeff,single(dT),zstart(1));


    for k=1:length(zstart)
        if photonratiofixed
           zstarth=repmat(zstart(k),1,size(imstack,3));
        else
            zstarth=zstart(k);
        end
         [Ph,CRLBh,LogLh]=fitterh(imstack,int32(5),shared,int32(iterations),splinecoeff,single(dT),zstarth);
        if first %first round
            P=Ph;
            CRLB=CRLBh;
            LogL=LogLh;
%             LogLr=ones(size(LogLh));
            first=false;
            color=zeros(size(LogLh))+phot;
            zpos=ones(size(LogLh));
        else
            indbetter=LogLh-LogL>0; %copy only everything if LogLh increases by more than rounding error.
            P(indbetter,:)=Ph(indbetter,:);
            CRLB(indbetter,:)=CRLBh(indbetter,:);
%             LogLr(indbetter)=LogL(indbetter)./LogLh(indbetter);
%             LogLr(~indbetter)=min(LogLr(~indbetter),LogL(~indbetter)./LogLh(~indbetter));
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
     [P,CRLB,LogL]=fitterh(imstack,int32(5),shared,int32(iterationsin),splinecoeff,single(dT),zstarth);
else
    color=[];
end

    ml=max(LogLp,[],2);
    LLsecond=quantile(ml./LogLp,1-1/length(PhotonRatios),2); 
%  clear(allfittersnames{fitter})
end


function  dT=ratiochannelshift(channelshift,PhotonRatio)
dS1 = [1, 1 ;1, 1 ;1, 1;1, PhotonRatio;1, 1];
npar=size(channelshift,1);
Nfits=size(channelshift,3);
noChannels=size(channelshift,2);
temp1 = repmat(dS1,[1 1 Nfits]);
dT = zeros(npar,noChannels*2,Nfits);
dT(:,1:noChannels,:)=channelshift;
dT(:,noChannels+1:2*noChannels,:)=temp1;
end
 
