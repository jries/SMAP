PSF=130;
pixel=130;
N=logspace(0,5,10);
Bg=[0 1 3 5 10 20 50 100 200];
figure(88);
% subplot(2,1,1)
% hold off
% subplot(2,1,2)
hold off
l={};
for k=1:length(Bg)
cmosn=1.3;
QE=.95;
lp1=Mortensen(QE*N,QE*Bg(k),PSF,pixel,cmosn);

cmosn=.8;
QE=.82;
lp2=Mortensen(QE*N,QE*Bg(k),PSF,pixel,cmosn);

% figure(89);
% subplot(2,1,1)
% loglog(N,lp1,'r',N,lp2,'b')
% hold on
% title('95 red, 82 blue')


vr=(lp2-lp1)./lp2;
% subplot(2,1,2)
semilogx(N,vr)
hold on

l{k}=num2str(Bg(k));
end
plot(N,0*N,'k')
legend(l)



function lp=Mortensen(N,Bg,PSF, pixel,cmosn)
b=sqrt(Bg+cmosn^2);

PSFa=sqrt(PSF^2+pixel^2/12);
v=PSFa^2./N.*(16/9+8*pi*PSFa^2*b^2./N/pixel^2);
lp=sqrt(v);
end
%mortensen