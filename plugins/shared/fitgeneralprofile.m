function [fitp,fitprof,fittext]=fitgeneralprofile(profile,x,p,sigma)
if ~isfield(p,'restrictsigma')
    p.restrictsigma=false;
end
if p.restrictsigma
    fif=@convoluteshape;
else
    fif=@convoluteshapefits;
end
switch p.fitmodel.Value
    case 1%Gauss
        [fitp,fitprof,fittext]=fitgauss(profile,x);
        fitp(5)=fitp(3);
    case 2%tophat
        [fitp,fitprof,fstart]=fif(x,profile,sigma,1);
        fittext=['Step L: ' 9 num2str(fitp(3),3) ];
    case 4%ring
        [fitp,fitprof,fstart]=fif(x,profile,sigma,3);
        fittext=['Ring R: ' 9 num2str(fitp(3),3) ];
    case 3%disk
        [fitp,fitprof,fstart]=fif(x,profile,sigma,2);
        fittext=['Disk R: ' 9 num2str(fitp(3),3)];   
    case 5 %double Gauss distance
%         [fitp,fitprof,fstart]=fif(x,profile,sigma,4);
        [fitp,fitprof,fittext]=fit2gauss(profile,x);
%         fittext=['Distance d: ' 9 num2str(fitp(4),3)];         

end
fittext={fittext};
if p.fitmodel.Value>1
if p.restrictsigma
    fittext(end+1)={['sigma: ' 9 num2str(sigma,4)]};
else
    fittext(end+1)={['sigma: ' 9 num2str(fitp(5),4)]}; 
end
end
    
function [fitp,fitprof,fittext]=fitgauss(profile,x)
[~,s]=getFWHM(profile,x);
s=s/2.6*(x(2)-x(1));
[mp, ip]=max(profile);
startp=[mp x(ip) s 0];
fitp=mygaussfit(x,profile,startp);
fitprof=mygaussforfit(fitp,x);
fittext=['sigma: ' 9 num2str(fitp(3),4)];

function [fitp,fitprof,fittext]=fit2gauss(profile,x)
%try with two Gaussian fits.
fp1=fit(x',profile','gauss1','Lower',[0 -inf 0]);
fp2=fit(x',profile'-(fp1(x)),'gauss1','Lower',[0 -inf 0]);

ft=fittype('a1*exp(-((x-b1)/c1)^2) + a2*exp(-((x-b2)/c1)^2)+d');
% fp=fit(x',profile','gauss2','Lower',[0 -inf 0 0 -inf 0]);
% fitp=[fp.a1 fp.b1 fp.c1 fp.a2 fp.b2 fp.c2];
startp=[fp1.a1 fp2.a1 fp1.b1 fp2.b1 fp2.c1 0];

fp=fit(x',profile',ft,'Lower',[0 0 -inf -inf 0 0],'StartPoint',startp);
% [~,s]=getFWHM(profile,x);
% s=s/2.6*(x(2)-x(1));
% [mp, ip]=max(profile);
% startp=[mp x(ip) s 0];
% fitp=mygaussfit(x,profile,startp);
fitp=[fp.a1 fp.a2 fp.b1 fp.b2-fp.b1 fp.c1 fp.d ];
fitprof=fp(x);
fittext=['Distance d: ' 9 num2str(abs(fp.b2-fp.b1),4)];

function [fwhm,fwhmind]=getFWHM(profile,x)
    [mp, ip]=max(profile);
    i1=find(profile(1:ip)>mp/2,1,'first');
    i2=find(profile(ip:end)>mp/2,1,'last')+ip-1;
    if isempty(i2)||isempty(i1)
        fwhm=[];
        fwhmind=1;
    else
        if i1==i2
            i1=i1-1;i2=i2+1;
        end
        
    fwhm=x(i2)-x(i1);
    fwhmind=i2-i1;
    end

