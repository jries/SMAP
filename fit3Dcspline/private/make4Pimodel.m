function [I,A,B,Ia,Aa,Ba]=make4Pimodel(allPSFso,phaseshifts,frequency,normf)
%re-weight every PSF by relative transmission?
for k=size(allPSFso,4):-1:1
    allPSFs(:,:,:,k)=allPSFso(:,:,:,k)/normf(k);
end

I1=(allPSFs(:,:,:,1)+allPSFs(:,:,:,3))/2;
I2=(allPSFs(:,:,:,2)+allPSFs(:,:,:,4))/2;
Iall=(I1+I2)/2;

z=(1:size(allPSFs,3))'-round(size(allPSFs,3)/2);
[A12,B12]=makeAB(allPSFs(:,:,:,1),allPSFs(:,:,:,2),Iall,z,frequency,phaseshifts(1),phaseshifts(2));
[A23,B23]=makeAB(allPSFs(:,:,:,2),allPSFs(:,:,:,3),Iall,z,frequency,phaseshifts(2),phaseshifts(3));
[A34,B34]=makeAB(allPSFs(:,:,:,3),allPSFs(:,:,:,4),Iall,z,frequency,phaseshifts(3),phaseshifts(4));
[A41,B41]=makeAB(allPSFs(:,:,:,4),allPSFs(:,:,:,1),Iall,z,frequency,phaseshifts(4),phaseshifts(1));
A=(A12+A23+A34+A41)/4;
B=(B12+B23+B34+B41)/4;
I=Iall;
if nargout>3 %also pass on individual calculations
    Aa(:,:,:,1)=A12;Aa(:,:,:,2)=A23;Aa(:,:,:,3)=A34;Aa(:,:,:,4)=A41;
    Ba(:,:,:,1)=B12;Ba(:,:,:,2)=B23;Ba(:,:,:,3)=B34;Ba(:,:,:,4)=B41;
    Ia(:,:,:,1)=I1;Ia(:,:,:,2)=I2;
end
end

function [A,B]=makeAB(P1,P2,I,z,frequency,phase1,phase2)
    A=zeros(size(I));B=zeros(size(I));
    for k=1:length(z)
        a1=2*frequency*z(k)+phase1;
        a2=2*frequency*z(k)+phase2;
        A(:,:,k)=(sin(a1).*(P2(:,:,k)-I(:,:,k))-sin(a2).*(P1(:,:,k)-I(:,:,k)))./(cos(a2).*sin(a1)-cos(a1).*sin(a2));
        B(:,:,k)=(-cos(a1).*(P2(:,:,k)-I(:,:,k))+cos(a2).*(P1(:,:,k)-I(:,:,k)))./(cos(a2).*sin(a1)-cos(a1).*sin(a2));
    end
end