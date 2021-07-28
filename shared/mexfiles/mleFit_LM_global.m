function [P,CRLB,LogL]=mleFit_LM_global(varargin)
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
P=[];CRLB=[];LogL=[]; %

%determine of it runs on GPU, otherwise use CPU as default
persistent fitter
allfitters={@GPUmleFit_LM_MultiChannel_noRestrict,@CPUmleFit_LM_MultiChannel};
allfittersnames={'GPUmleFit_LM_MultiChannel_noRestrict','CPUmleFit_LM_MultiChannel'};
if isempty(fitter)
    for k=1:length(allfitters)
        try
%             allfitters{k}(ones(7,7,1,2,'single'),single([1 1 1 0 0]),varargin{4},zeros(5,2));
            allfitters{k}(varargin{1:5});
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
imstack=single(varargin{1});
shared=int32(varargin{2});
splinecoeff=single(varargin{4});
channelshift=single(varargin{5});

coeffsize=size(splinecoeff);
fitmode=5;
% if nargin<7||isempty(varargin{7})
%     varmap=0; %emccd
% end 

if nargin<7 || isempty(varargin{7})
    PhotonRatios=[];
else
    PhotonRatios=varargin{7};
    shared(4,:)=1; %if we fit with ratios, N needs to be linked
end

if nargin<6||isempty(varargin{6})
    zstart=single(coeffsize(3)/2); %emccd
else
    zstart=single(coeffsize(3)/2+varargin{6});
end



if nargin>2 && ~isempty(varargin{3})
    iterations=varargin{3};
else
    iterations=30;
end


%Yiming: PhotonRatios is a vector with photon ratios. You know best how to
%pass it on to the function. If length(PhotonRatios)>1: loop over, choose
%highest LL, save color (e.g. at P(:, end+1))


[P,CRLB,LogL]=allfitters{fitter}(imstack,shared,iterations,splinecoeff,channelshift,zstart(1));

if length(zstart)>1
    for k=2:length(zstart)
        [Ph,CRLBh,LogLh]=allfitters{fitter}(imstack,shared,iterations,splinecoeff,channelshift,zstart(k));
%         indbettero=LogLh<LogL;
        indbetter=LogLh-LogL>1e-4; %copy only everything if LogLh increases by more than rounding error.
        P(indbetter,:)=Ph(indbetter,:);
        CRLB(indbetter,:)=CRLBh(indbetter,:);
        LogL(indbetter)=LogLh(indbetter);
    end
end


 clear(allfittersnames{fitter})

