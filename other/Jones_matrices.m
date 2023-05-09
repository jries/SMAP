%Jones martrices
v=[1;0];
vstart=linear(pi/2);
clear phases dp
phases(1,1,:)=-pi/2:0.01:pi/2;


[vout, dangle]=EomSLM(vstart,phases,0);
[vout1, dangle1]=EomSLM(vstart,phases,pi/2*1.9);
[vout2, dangle2]=EomSLM(vstart,phases,-pi/2*1.9);
sphases=squeeze(phases);
outchannel=1;
figure(88);subplot(1,2,1);hold off;

plot(sphases,squeeze(dangle),'r',sphases,squeeze(abs(vout(outchannel,:,:)).^2/2),'b',sphases,0*sphases,'k')
legend('phase diff','intensity(0)')
xlabel('phase shift')
hold on
plot(sphases,squeeze(dangle1),'r--',sphases,squeeze(abs(vout1(outchannel,:)).^2/2),'b--')
plot(sphases,squeeze(dangle2),'r--',sphases,squeeze(abs(vout2(outchannel,:)).^2/2),'b--')
title('EOM delay for 3 lateral phases')

dp(1,1,:)=-2*pi:0.1:2*pi;
[voutd, dangled]=EomSLM(vstart,-pi/4,dp);
[voutd1, dangled1]=EomSLM(vstart,0,dp);
[voutd2, dangled2]=EomSLM(vstart,pi/4,dp);
sphases=squeeze(dp);

figure(88);subplot(1,2,2);hold off;

plot(sphases,squeeze(dangled),'r',sphases,squeeze(abs(voutd(outchannel,:,:)).^2/2),'b',sphases,0*sphases,'k')
legend('phase diff','intensity(0)')
xlabel('lateral position')
hold on
plot(sphases,squeeze(dangled1),'r--',sphases,squeeze(abs(voutd1(outchannel,:)).^2/2),'b--')
plot(sphases,squeeze(dangled2),'r--',sphases,squeeze(abs(voutd2(outchannel,:)).^2/2),'b--')

title('lateral phases for 3 EOM delays')

function [vout, dangle]=EomSLM(vstart,phases,dphi)
if nargin<3
    dphi=0;
end
vphase=pagemtimes(phaseplate(phases),vstart);
 % vphase=pagemtimes(QWP(phases),vstart);
vphase=pagemtimes(addphase(dphi),vphase);
vrot=pagemtimes(addphase(-dphi),pagemtimes(HWP45,vphase));
vsum=(vrot+vphase);
vout=pagemtimes(polarizer(0),vsum);
outchannel=1;
dangle=mod(angle(vrot(outchannel,:))-angle(vphase(outchannel,:))+pi,2*pi)-pi;
end

function out=waveplate(phase,angle)
out=exp(-1i*phase/2).*...
    [cos(angle).^2+exp(1i*phase)*sin(angle).^2 (1-exp(1i*phase))*cos(angle).*sin(angle)
    (1-exp(1i*phase))* cos(angle).*sin(angle) sin(angle).^2+exp(1i*phase)*cos(angle).^2];
end

function out=QWP(angle)
out=waveplate(pi/2,angle);
end
function out=HWP(angle)
out=waveplate(pi,angle);
end
function out=HWP45
out=HWP(pi/4);
end
function out=QWP45
out=QWP(pi/4);
end
function out=phaseplate(phase)
out=waveplate(2*phase,0);
end
function out=addphase(phase)
out=exp(1i*phase);
end

function out=linear(angle)
v=[1;0];
out=1i*HWP(angle/4)*v;
end

function out=polarizer(angle)
out=[cos(angle)^2 cos(angle)*sin(angle)
   cos(angle)*sin(angle) sin(angle)^2];
end