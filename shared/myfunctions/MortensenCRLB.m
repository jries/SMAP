function [lp,errphot]=MortensenCRLB(N,Bg,PSF, pixel,cmosn)
%N, Bg in photons, cmosn in e- (std, not variance). 
b=sqrt(Bg.^2+cmosn^2);

PSFa=sqrt(PSF.^2+pixel^2/12);
v=PSFa.^2./N.*(16/9+8*pi*PSFa.^2.*b./N/pixel^2);  %here b is background photons per pixel, in Mortensen it is b^2, as b is the noise
lp=sqrt(v/2); %the Formula takes into account EMCCD and doubling ouf noise because of that.

s_a=PSF/pixel; %sigmapsf/pixelsize
tau=2*pi*(b).*(s_a.^2+1/12)./N;
errphot2=N.*(1+4*tau+sqrt(tau./(14*(1+2*tau)))); %This is Rieger...
errphot=sqrt(errphot2);
end



% function [lp,errphot]=Mortensen(N,Bg,PSF, pixel,cmosn)
% b=sqrt(Bg+cmosn^2);
% 
% PSFa=sqrt(PSF^2+pixel^2/12);
% v=PSFa^2./N.*(16/9+8*pi*PSFa^2*b.^2./N/pixel^2);
% lp=sqrt(v);
% 
% 
% s_a=PSF/pixel; %sigmapsf/pixelsize
% tau=2*pi*(Bg)*(s_a^2+1/12)./N;
% errphot2=N.*(1+4*tau+sqrt(tau./(14*(1+2*tau)))); %This is Rieger...
% errphot=sqrt(errphot2);
% end