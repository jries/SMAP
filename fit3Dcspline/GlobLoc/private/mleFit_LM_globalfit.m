function [P,CRLB,LogL]=mleFit_LM_globalfit(varargin)
% varargin:
%imstack, sharedflag, iterations, spline coefficients, channelshift,
%fitmode, varmap,zstart

% 1. imagestack (single)
% 2. fittype: 1 for 2D Gauss, 2 for spline
% 3. shared parameters

% 4. iterations (default=50)
% 5. paramters for fitters:

%   cspline: cspline coefficients (single)
%   Gauss :init sigma
% 6. shifts dT
% 7. scmoS varmap
% 8. silent (suppress output if 1)

% 9. z start parameter (more than one: return solution with maximum
%       LIkelihood). Units: distance of stack calibration, center based


%Output:
%P
%1. X, Y, Photons, Background, PSFxy, Iterations
%2. X, Y, Photons, Z, Background, PSFxy, Iterations

%CRLB: cramer-rao lower bounds, as in P
%LogL: log-likelihood.

%Only for fitmode 6: P1 etc: results with z-startparameter<0, P2 etc:
%results with z-startparameter>0

% [P,CRLB, LL] =CPUmleFit_LM_MultiChannel(fitstack,int32(sharedA),50,coeffh,single(dT),fitmode);
P=[];CRLB=[];LogL=[]; %

if nargin<6
    disp('At least 6 input parameters');
    return
end
imstack=single(varargin{1});
fittype = int32(varargin{2});
shared=int32(varargin{3});
iterations=int32(varargin{4});
splinecoeff=single(varargin{5});
channelshift=single(varargin{6});

coeffsize=size(splinecoeff);
Nfits = size(imstack,3);

%determine of it runs on GPU, otherwise use CPU as default
persistent fitter
allfitters={@GPUmleFit_LM_MultiChannel,@CPUmleFit_LM_MultiChannel};
allfittersnames={'GPUmleFit_LM_MultiChannel','CPUmleFit_LM_MultiChannel'};
if isempty(fitter)
    for k=1:length(allfitters)
        try
%             allfitters{k}(ones(7,7,1,2,'single'),single([1 1 1 0 0]),varargin{4},zeros(5,2));
            allfitters{k}(imstack,fittype,shared,iterations,splinecoeff,channelshift);
            fitter=k;
            break
        catch err
            % fitter did not work
        end
    end
    disp(['using: ' char(allfitters{fitter})]);
end

if nargin>6 && ~isempty(varargin{7})
    varmap=single(varargin{7});
else
    varmap=0;
end

if nargin>7 && ~isempty(varargin{8})
    silent=single(varargin{8});
else
    silent=1;
end

if nargin>8 && ~isempty(varargin{9})&& (fittype==2)
    zstart=single(coeffsize(3)/2+varargin{9});
elseif (fittype==2)
    zstart=repmat(single(coeffsize(3)/2),[Nfits,1]);
else
    zstart=0;
end


[P,CRLB,LogL]=allfitters{fitter}(imstack,fittype,shared,iterations,splinecoeff,channelshift,varmap,silent,zstart(:,1));


if size(zstart,2)>1
    for k=2:size(zstart,2)
        [Ph,CRLBh,LogLh]=allfitters{fitter}(imstack,fittype,shared,iterations,splinecoeff,channelshift,varmap,silent,zstart(:,k));
%         indbettero=LogLh<LogL;
        indbetter=LogLh-LogL>1e-4; %copy only everything if LogLh increases by more than rounding error.
        P(indbetter,:)=Ph(indbetter,:);
        CRLB(indbetter,:)=CRLBh(indbetter,:);
        LogL(indbetter)=LogLh(indbetter);
    end
end
%     [P,CRLB,LogL]=allfitters{fitter}(imstack,shared,iterations,splinecoeff,channelshift,zstart);
% end
%%
% 
% <latex>
% \begin{tabular}{|c|c|} \hline
% $n$ & $n!$ \\ \hline
% 1 & 1 \\
% 2 & 2 \\
% 3 & 6 \\ \hline
% \end{tabular}
% </latex>
% 
 clear(allfittersnames{fitter})

