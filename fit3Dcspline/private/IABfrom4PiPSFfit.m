function out=IABfrom4PiPSFfit(PSF, phaseshift,frequency,roisizexy,numframes,zshift0)
if nargin<6
    zshift0 =[0 0 0 0];
end
% startpar=[dx dy dz phase norm] all vectors length 3 ;
% dx=par(1:3);
% dy=par(4:6);
% dz=par(7:9);
% normf=par(10:12);
%shift PSF by dx,dy,dz

s=size(PSF);
rx=ceil((roisizexy+1)/2)+2;mp=ceil((s(1)+1)/2);mpz=ceil((s(3)+1)/2);rz=min(numframes,mpz-2);
PSFs=PSF(mp-rx:mp+rx,mp-rx:mp+rx,mpz-rz:mpz+rz,:);

% startpar=[0 0 0 0 0 0 0 0 0 1 1 1];
startpar=[0 0 0 0 0 0  1 1 1];
% startpar=[0 0 0 0 0 0  1 1 1 phaseshift]
% startpar=[0 0 0 0 0 0]
% fixpar=[phaseshift frequency];
zshift0h=zshift0(2:4)-zshift0(1);
fixpar=[phaseshift frequency zshift0h];

[PSFstart,mstart]=recoverPSF(startpar,PSFs,fixpar);

options=optimoptions('lsqnonlin','Display','iter');
fitpar=lsqnonlin(@errPSFIAB,startpar,[],[],options,PSF,fixpar);
% options=optimset('lsqcurvefit');
% options.FinDiffRelStep=.01;
% PSFc=PSFs(3:end-2,3:end-2,3:end-2,:);


% fitpar=lsqcurvefit(@recoverPSF,startpar,PSFs,PSFc,[],[],options,fixpar);

[PSFrecovered,out]=recoverPSF(fitpar,PSF,fixpar);

out.dx=[0 fitpar(1:3)];
out.dy=[0 fitpar(4:6)];
out.dz=[0 zshift0h];
% out.dz=[0 fitpar(7:9)];
% out.normf=[1 fitpar(10:12)];
if length(fitpar)>6
out.normf=[1 fitpar(7:9)];
else
    out.normf=ones(4,1);
end
if length(fitpar)>9
    phaseshift=fitpar(10);
    if phaseshift>pi
        phaseshift=phaseshift-2*pi;
    end
end
out.frequency=frequency;
out.phaseshifts=[-pi phaseshift 0 phaseshift+pi];

end

function err=errPSFIAB(par,PSF,fixpar)
[PSFmSs,out,dPSF,dI,dA,dB]=recoverPSF(par,PSF,fixpar);
err=vertcat(dPSF(:),dI(:),dA(:),dB(:));
end

function [PSFmSs,out,dPSF,dI,dA,dB]=recoverPSF(par,PSF,fixpar)

dx=par(1:3);
dy=par(4:6);
% dz=par(7:9);
% normf=par(10:12);
if length(par)>6
normf=par(7:9);
else
    normf=ones(1,3);
end
if length(par)>9
    phaseshift=par(10);
else
phaseshift=fixpar(1);
end
frequency=fixpar(2);
dz=fixpar(3:5);

phaseshifts=[-pi phaseshift 0 phaseshift+pi];

%shift PSF and rescale
xn=1:size(PSF,1);yn=1:size(PSF,2);zn=1:size(PSF,3);
[Xq,Yq,Zq]=meshgrid(yn,xn,zn);
PSFS=zeros(size(PSF));
PSFS(:,:,:,1)=PSF(:,:,:,1);
for k=2:4
    PSFS(:,:,:,k)=interp3(PSF(:,:,:,k),Xq-dx(k-1),Yq-dy(k-1),Zq-dz(k-1),'cubic',0);%/normf(k-1);
end

% calculate IAB
[I,A,B,Ii,Ai,Bi]=make4Pimodel(PSFS,phaseshifts,frequency,[1 normf]);

s=size(PSF);
rim=2;
r1=rim+1:s(1)-rim;r2=rim+1:s(2)-rim;r3=rim+1:s(3)-rim;
dA=zeros(s(1)-2*rim,s(2)-2*rim,s(3)-2*rim,4);dI=zeros(s(1)-2*rim,s(2)-2*rim,s(3)-2*rim,2);dB=dA;
dPSF=zeros(s(1)-2*rim,s(2)-2*rim,s(3)-2*rim,4);

for k=1:4
    dA(:,:,:,k)=A(r1,r2,r3)-Ai(r1,r2,r3,k);
    dB(:,:,:,k)=B(r1,r2,r3)-Bi(r1,r2,r3,k);
end
for k=1:2
    dI(:,:,:,k)=I(r1,r2,r3)-Ii(r1,r2,r3,k);
end


%make PSF again
PSFm=makePSF(I,A,B,frequency, phaseshifts, [1 normf]);

%shift back for fitting
PSFmS=zeros(size(PSF));
PSFmS(:,:,:,1)=PSFm(:,:,:,1);
for k=2:4
    PSFmS(:,:,:,k)=interp3(PSFm(:,:,:,k),Xq+dx(k-1),Yq+dy(k-1),Zq+dz(k-1),'cubic',0);
end

for c=1:4
    dPSF(:,:,:,c)=(PSFmS(r1,r2,r3,c)-PSF(r1,r2,r3,c));
end

PSFmSs=PSFmS(3:end-2,3:end-2,3:end-2,:);
% PSFm=PSFm(2:end-1,2:end-1,2:end-1,2:end-1);
out.PSF=PSFmS;
out.I=I;
out.A=A;
out.B=B;
end

% 
% function [I,A,B]=make4Pimodel(allPSFso,phaseshifts,frequency,normf)
% %re-weight every PSF by relative transmission?
% for k=size(allPSFso,4):-1:1
%     allPSFs(:,:,:,k)=allPSFso(:,:,:,k)/normf(k);
% end
% 
% I1=(allPSFs(:,:,:,1)+allPSFs(:,:,:,3))/2;
% I2=(allPSFs(:,:,:,2)+allPSFs(:,:,:,4))/2;
% Iall=(I1+I2)/2;
% 
% z=(1:size(allPSFs,3))'-round(size(allPSFs,3)/2);
% [A12,B12]=makeAB(allPSFs(:,:,:,1),allPSFs(:,:,:,2),Iall,z,frequency,phaseshifts(1),phaseshifts(2));
% [A41,B41]=makeAB(allPSFs(:,:,:,4),allPSFs(:,:,:,1),Iall,z,frequency,phaseshifts(4),phaseshifts(1));
% [A23,B23]=makeAB(allPSFs(:,:,:,2),allPSFs(:,:,:,3),Iall,z,frequency,phaseshifts(2),phaseshifts(3));
% [A34,B34]=makeAB(allPSFs(:,:,:,3),allPSFs(:,:,:,4),Iall,z,frequency,phaseshifts(3),phaseshifts(4));
% A=(A12+A23+A34+A41)/4;
% B=(B12+B23+B34+B41)/4;
% I=Iall;
% end

function [A,B]=makeAB(P1,P2,I,z,frequency,phase1,phase2)
    A=zeros(size(I));B=zeros(size(I));
    for k=1:length(z)
        a1=2*frequency*z(k)+phase1;
        a2=2*frequency*z(k)+phase2;
        A(:,:,k)=(sin(a1).*(P2(:,:,k)-I(:,:,k))-sin(a2).*(P1(:,:,k)-I(:,:,k)))./(cos(a2).*sin(a1)-cos(a1).*sin(a2));
        B(:,:,k)=(-cos(a1).*(P2(:,:,k)-I(:,:,k))+cos(a2).*(P1(:,:,k)-I(:,:,k)))./(cos(a2).*sin(a1)-cos(a1).*sin(a2));
    end
end

function PSFo=makePSF(I,A,B,frequency, phaseshifts, normf)
s=size(I);
PSF=zeros(s(1)*s(2),s(3),4);
z=(1:s(3))-round(s(3)/2);
Ir=reshape(I,s(1)*s(2),s(3));
Br=reshape(B,s(1)*s(2),s(3));
Ar=reshape(A,s(1)*s(2),s(3));
for k=1:4
    PSF(:,:,k)=normf(k)*(Ir+Ar.*cos(2*frequency*z+phaseshifts(k))+Br.*sin(2*frequency*z+phaseshifts(k)));
end

PSFo=reshape(PSF,s(1),s(2),s(3),4);
end