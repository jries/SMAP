function [P,CRLB,LogL]=mleFit_LM_global_gauss(varargin)
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
% 6. silent (suppress output if 1)
% 7. scmoS varmap
% 8. z start parameter (more than one: return solution with maximum
%       LIkelihood). Units: distance of stack calibration, center based


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
allfitters={@GPUmleFit_LM_MultiChannel_Gauss,@CPUmleFit_global_Gauss};
allfittersnames={'GPUmleFit_LM_MultiChannel_Gauss','CPUmleFit_global_Gauss'};

% allfitters={@CPUmleFit_global_Gauss};
% allfittersnames={'CPUmleFit_global_Gauss'};
if isempty(fitter)
    for k=1:length(allfitters)
        try
%             allfitters{k}(ones(7,7,1,2,'single'),single([1 1 1 0 0]),varargin{4},zeros(5,2));
            allfitters{k}(varargin{:});
            fitter=k;
            break
        catch err
            % fitter did not work
        end
    end
    disp(['using: ' char(allfitters{fitter})]);
end

[P,CRLB,LogL]=allfitters{fitter}(varargin{:});

 clear(allfittersnames{fitter})

