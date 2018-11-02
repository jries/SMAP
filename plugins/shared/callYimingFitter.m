function [P,CRLB,LogL]=callYimingFitter(varargin)
sdafjd
%imstack, startpsf/coeff, iterations, fitmode, isemccd, hidereport
persistent fitter
% cpufitter=@kernel_MLEfit_Spline_LM_SMAP_v2_nointerp;
% cpufitter=@CPUmleFit_LM;
% gpufitter=@GPUmleFit_LM;
allfitters={@GPUmleFit_LM,@GPUmleFit_LM_CUDA75,@CPUmleFit_LM,@kernel_MLEfit_Spline_LM_SMAP_v2_nointerp};
allfittersnames={'GPUmleFit_LM','GPUmleFit_LM_CUDA75','CPUmleFit_LM','kernel_MLEfit_Spline_LM_SMAP_v2_nointerp'};
if isempty(fitter)
    for k=1:length(allfitters)
        try
            allfitters{k}(ones(7,'single'),1,10,1,0,0);
            fitter=k;
            break
        end
    end
    disp(['using: ' char(allfitters{fitter})]);
%     try
%         gpufitter(ones(7,'single'),1,10,1,0);
%         fitter=2;
%     catch
%         disp('GPU didnt work, run on CPU (slow)');
%         fitter=1;
%     end
end
varargingpu=varargin;
if length(varargingpu)==5
    varargingpu{6}=1;
end
% switch fitter
%     case 1
%         [P]=cpufitter(varargin{:});
%     case 2
% 
%         [P,CRLB,LogL]=gpufitter(varargingpu{:});
% end

[P,CRLB,LogL]=allfitters{fitter}(varargingpu{:});
 clear(allfittersnames{fitter})

%
