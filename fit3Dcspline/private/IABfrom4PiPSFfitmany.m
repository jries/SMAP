function out=IABfrom4PiPSFfitmany(img,startp, phaseshift,frequency,roisizexy,numframes)
% startpar=[dx dy dz phase norm] all vectors length 3 ;
% par= dxi
%      dyi
%      dzi
%      phii
%      Ni
%       %or rather normp for the four quadrants?
%      dxp
%      dyp
%      normp

% normf=par(10:12);
%shift PSF by dx,dy,dz
PSF=img.imstack;

s=size(PSF);
rx=ceil((roisizexy+1)/2)+2;mp=ceil((s(1)+1)/2);mpz=ceil((s(3)+1)/2);rz=min(numframes,mpz-1);
PSFs=PSF(mp-rx:mp+rx,mp-rx:mp+rx,mpz-rz:mpz+rz,:,:);


Nistart=startp.N;
% dxstart=zeros(1,s(4));
% dxstart=startp.x0-startp.x0(1);
% dystart=startp.y0-startp.y0(1);
PSFg=img.fit.PSF;
dxstart=startp.x0-startp.x0(1)+PSFg.dx;
dystart=startp.y0-startp.y0(1)+PSFg.dy;
dzstart=startp.z0-startp.z0(1);
phasestart=startp.phase-startp.phase(1);


% dxpstart=-PSFg.dx;
% dypstart=-PSFg.dy;
% normpstart=PSFg.normf;
normpstart=[1 1 1 1]; 

% startpar=[-dxstart(:);-dystart(:);-dzstart(:);phasestart(:); Nistart(:); -dxpstart(:);-dypstart(:);normpstart(:)];
% startpar=[-dxstart(:);-dystart(:);-dzstart(:);phasestart(:); Nistart(:); -dxpstart(:);-dypstart(:);normpstart(:); phaseshift];
startpar=[-dxstart(:);-dystart(:);-dzstart(:);phasestart(:); Nistart(:);normpstart(:); phaseshift];
% startpar=[dxstart(:);dxstart(:);dxstart(:);dxstart(:); Nistart(:); dxpstart(:);dypstart(:);normpstart(:);phaseshift;frequency];
% startpar=[0 0 0 0 0 0 0 0 0 1 1 1];
% startpar=[0 0 0 0 0 0  1 1 1];
% fixpar=[phaseshift frequency];
% zshift0h=zshift0(2:4)-zshift0(1);
fixpar=[phaseshift frequency];
err=errPSFIAB(startpar,PSFs,fixpar);

sum(err.^2)

% options=optimoptions('fminunc','Display','iter');
% fitpar=fminunc(@errPSFIAB,startpar,options,PSF,fixpar);

options=optimoptions('lsqnonlin','Display','iter');
fitpar=lsqnonlin(@errPSFIAB,startpar,[],[],options,PSF,fixpar);

[PSFstart,mstart]=recoverPSF(startpar,PSFs,fixpar);

options=optimoptions('lsqcurvefit','Display','iter');
options.FinDiffRelStep=.01;
PSFc=PSFs(3:end-2,3:end-2,3:end-2,:,:);


fitpar=lsqcurvefit(@recoverPSF,startpar,PSFs,PSFc,[],[],options,fixpar);

[PSFrecovered,out]=recoverPSF(fitpar,PSF,fixpar);

out.dx=[0 fitpar(1:3)];
out.dy=[0 fitpar(4:6)];
out.dz=[0 zshift0h];
% out.dz=[0 fitpar(7:9)];
% out.normf=[1 fitpar(10:12)];
out.normf=[1 fitpar(7:9)];
out.frequency=frequency;
out.phaseshifts=[-pi phaseshift 0 phaseshift+pi];

end

function err=errPSFIAB(par,PSF,fixpar)
[PSFmSs,out,dPSF,dI,dA,dB]=recoverPSF(par,PSF,fixpar);
err=vertcat(dPSF(:),dI(:),dA(:),dB(:));


%later: maybe calculate a bead-dependent weighting factor for 'robust'
%fitting.

% erro=sum(err.^2);
%calculate an error term for the optimizer to minimize
% final model I, A, B, PSF
% for each bead 4x A,B, 2x I, 1x PSF: just join?
%call recover PSF, also pass out all IAB

end

function [PSFmSs,out,dPSF,dI,dA,dB]=recoverPSF(par,PSF,fixpar)
rim=4;
s=size(PSF);
numbeads=s(4);

% dx=par(1:numbeads);
% dy=par(numbeads+1:2*numbeads);
% dz=par(2*numbeads+1:3*numbeads);
% phi=par(3*numbeads+1:4*numbeads);
% Ni=par(4*numbeads+1:5*numbeads);


dx=reshape(par(1:numbeads*4),numbeads,4);
dy=reshape(par(numbeads*4+1:8*numbeads),numbeads,4);
dz=par(8*numbeads+1:9*numbeads);
phi=par(9*numbeads+1:10*numbeads);
Ni=par(10*numbeads+1:11*numbeads);

% off=5*numbeads;
off=11*numbeads;
% dxp=par(off+1:off+4);
% dyp=par(off+5:off+8);
normp=par(off+1:off+4);
phaseshift=par(off+5);
% frequency=par(off+14);

dx(1)=0;
dy(1)=0;
dz(1)=0;
phi(1)=0;
dxp(1)=0;
dyp(1)=0;
normp(1)=1;

% phaseshift=fixpar(1);
frequency=fixpar(2);

phaseshifts=[-pi phaseshift 0 phaseshift+pi];

%shift PSF and rescale
xn=1:size(PSF,1);yn=1:size(PSF,2);zn=1:size(PSF,3);
[Xq,Yq,Zq]=meshgrid(yn,xn,zn);
PSFS=zeros(size(PSF));
Ia=zeros(s(1),s(2),s(3),s(4));Aa=Ia;Ba=Ia;

dA=zeros(s(1)-2*rim,s(2)-2*rim,s(3)-2*rim,s(4),4);dI=zeros(s(1)-2*rim,s(2)-2*rim,s(3)-2*rim,s(4),2);dB=dA;
r1=rim+1:s(1)-rim;r2=rim+1:s(2)-rim;r3=rim+1:s(3)-rim;
dPSF=zeros(s(1)-2*rim,s(2)-2*rim,s(3)-2*rim,s(4),4);

for b=1:numbeads
    for c=1:4
        dxh=dx(b,c);
        dyh=dy(b,c);
%         dxh=dx(b)+dxp(c);
%         dyh=dy(b)+dyp(c);
        dzh=dz(b);
        PSFS(:,:,:,b,c)=interp3(PSF(:,:,:,b,c),Xq-dxh,Yq-dyh,Zq-dzh,'cubic',0);%/normf(k-1);
    end
        [I,A,B,Ii,Ai,Bi]=make4Pimodel(squeeze(PSFS(:,:,:,b,:)),phaseshifts+phi(b),frequency,normp*Ni(b));
        Ia(:,:,:,b)=I;
        Aa(:,:,:,b)=A;
        Ba(:,:,:,b)=B;
        for k=1:2
            dI(:,:,:,b,k)=Ia(r1,r2,r3)-Ii(r1,r2,r3,k);
        end
        for k=1:4
            dA(:,:,:,b,k)=Aa(r1,r2,r3)-Ai(r1,r2,r3,k);
            dB(:,:,:,b,k)=Ba(r1,r2,r3)-Bi(r1,r2,r3,k);
        end
%     end
end

% calculate IAB
% [I,A,B]=make4Pimodel(PSFS,phaseshifts,frequency,[1 normf]);

%make PSF again
I=mean(Ia,4);A=mean(Aa,4);B=mean(Ba,4);
PSFm=makePSF(I,A,B,frequency, phaseshifts,  normp);

%shift back for fitting
PSFmS=zeros(size(PSF));
% PSFmS(:,:,:,1)=PSFm(:,:,:,1);
% for k=2:4
%     PSFmS(:,:,:,k)=interp3(PSFm(:,:,:,k),Xq+dx(k-1),Yq+dy(k-1),Zq+dz(k-1),'cubic',0);
% end

for b=1:numbeads
    for c=1:4
        dxh=dx(b,c);
        dyh=dy(b,c);
%         dxh=dx(b)+dxp(c);
%         dyh=dy(b)+dyp(c);
        dzh=dz(b);
        PSFmS(:,:,:,b,c)=interp3(PSFm(:,:,:,c),Xq+dxh,Yq+dyh,Zq+dzh,'cubic',0)*Ni(b);%/normf(k-1);
        dPSF(:,:,:,b,c)=(PSFmS(r1,r2,r3,b,c)-PSF(r1,r2,r3,b,c))/Ni(b);
    end
end

PSFmSs=PSFmS(1+rim:end-rim,1+rim:end-rim,1+rim:end-rim,:,:);

% PSFm=PSFm(2:end-1,2:end-1,2:end-1,2:end-1);
out.PSF=PSFmS;
out.I=I;
out.A=A;
out.B=B;
out.dA=dA;out.dB=dB;out.dI=dI;out.dPSF=dPSF;
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