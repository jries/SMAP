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

%Only for fitmode 6: P1 etc: results with z-startparameter<0, P2 etc:
%results with z-startparameter>0

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

% if nargin<7||isempty(varargin{7})
%     varmap=0; %emccd
% end 
if nargin<9||isempty(varargin{9})
    zstart=single(coeffsize(3)/2); 
else
    zstart=single(coeffsize(3)/2+varargin{9});
end

phase=single(varargin{8});
if isempty(varargin{3})
    iterations=varargin{3};
else
    iterations=30;
end

% if fitmode==6
%     fitmode=5;
%     zstart=single([-coeffsize(3)/4, coeffsize(3)/4]+coeffsize(3)/2);
% %     zstart=single([-coeffsize(3)/3, coeffsize(3)/3]+coeffsize(3)/2);
% end
imstack=single(varargin{1});
shared=uint32(varargin{2});
channelshift=single(varargin{7});

% if fitmode==6 %2D fit: find proper results
% %     [P1,CRLB1,LogL1,P2,CRLB2,LogL2]=allfitters{fitter}(varargin{:});
% %     ind1=LogL1>=LogL2;
% %     ind2=LogL1<LogL2;
% %     P=zeros(size(P1),'single');CRLB=zeros(size(CRLB1),'single');LogL=zeros(size(LogL1),'single');
% %     P(ind1,:)=P1(ind1,:);P(ind2,:)=P2(ind2,:);
% %     CRLB(ind1,:)=CRLB1(ind1,:);CRLB(ind2,:)=CRLB2(ind2,:);
% %     LogL(ind1,:)=LogL1(ind1,:);LogL(ind2,:)=LogL2(ind2,:);
[P,CRLB,LogL]=allfitters{fitter}(imstack,shared,iterations,I,A,B,channelshift,phase,zstart(1));


if length(zstart)>1
    for k=2:length(zstart)
        [Ph,CRLBh,LogLh]=allfitters{fitter}(imstack,shared,iterations,I,A,B,channelshift,phase,zstart(k));
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

