%Jones martrices
v=[1;0];
vstart=linear(pi/2);

phases=-pi/2:0.01:pi/2;
% phases=[0 pi/4 pi];
clear vphase
for k=length(phases):-1:1
  vphase(:,k)=phaseplate(phases(k))*vstart;
end
vrot=HWP45*vphase;
vsum=(vrot+vphase);
vout=polarizer(0)*vsum;
outchannel=1;
dangle=mod(angle(vrot(outchannel,:))-angle(vphase(outchannel,:))+pi,2*pi)-pi;

figure(88);plot(phases,dangle,phases,abs(vout(outchannel,:)).^2/2,phases,0*phases)
legend('phase diff','intensity(0)')
xlabel('phase shift')



function out=waveplate(phase,angle)
out=exp(-1i*phase/2)*...
    [cos(angle)^2+exp(1i*phase)*sin(angle)^2 (1-exp(1i*phase))*cos(angle)*sin(angle)
    (1-exp(1i*phase))* cos(angle)*sin(angle) sin(angle)^2+exp(1i*phase)*cos(angle)^2];
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