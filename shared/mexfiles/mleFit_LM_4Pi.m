function [P,CRLB,LogL]=mleFit_LM_4Pi(varargin)
% varargin:
%imstack, sharedflag, iterations, spline coefficients, channelshift,
%fitmode, varmap,zstart

% 1. imagestack (single)
% 2. shared parameters

% optional:
% 3. iterations (default=50)
% 4. cspline I
% 5. cspline A
% 6. cspline B
% 7. shifts dT
% 8. phi0
% 9. z start parameter (more than one: return solution with maximum
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



% [P,CRLB, LL] =CPUmleFit_LM_MultiChannel(fitstack,int32(sharedA),50,coeffh,single(dT),fitmode);
% P=[];CRLB=[];LogL=[]; %
% 
% %determine of it runs on GPU, otherwise use CPU as default
% persistent fitter
allfitters={@CPUmleFit_LM_MultiChannel_4pi};
allfittersnames={'CPUmleFit_LM_MultiChannel_4pi'};
% if isempty(fitter)
%     for k=1:length(allfitters)
%         try
% %             allfitters{k}(ones(7,7,1,2,'single'),single([1 1 1 0 0]),varargin{4},zeros(5,2));
%             allfitters{k}(varargin{1:5});
%             fitter=k;
%             break
%         catch err
%             % fitter did not work
%         end
%     end
%     disp(['using: ' char(allfitters{fitter})]);
% end

fitter=1;

I=single(varargin{4});
A=single(varargin{5});
B=single(varargin{6});

coeffsize=size(I);

if nargin<10||isempty(varargin{10})
    p0=single(0); 
else
    p0=single(varargin{10});
end


if nargin<9||isempty(varargin{9})
    zstart=single(coeffsize(3)/2); 
else
    zstart=single(coeffsize(3)/2+varargin{9});
end

phase=single(varargin{8});
if ~isempty(varargin{3})
    iterations=varargin{3};
else
    iterations=30;
end

imstack=single(varargin{1});
shared=uint32(varargin{2});
channelshift=single(varargin{7});
if numel(zstart)~=size(imstack,3)
    zstart=ones(size(imstack,3),1,'single')*zstart;
end
if numel(p0)==1
    p0=p0*ones(size(imstack,3),1,'single');
end

[P,CRLB,LogL]=allfitters{fitter}(imstack,shared,iterations,I,A,B,channelshift,phase,zstart(:,1),p0(:));

sz0=size(zstart);
if sz0(2)>1
    for k=2:sz0(2)
        [Ph,CRLBh,LogLh]=allfitters{fitter}(imstack,shared,iterations,I,A,B,channelshift,phase,zstart(:,k),p0(:));
%         indbettero=LogLh<LogL;
        indbetter=LogLh-LogL>1e-4; %copy only everything if LogLh increases by more than rounding error.
        P(indbetter,:)=Ph(indbetter,:);
        CRLB(indbetter,:)=CRLBh(indbetter,:);
        LogL(indbetter)=LogLh(indbetter);
    end
end
 clear(allfittersnames{fitter})

