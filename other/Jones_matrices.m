%Jones martrices
v=[0;1];
 vlin45=linear(-pi/4);  


vstart=v;
clear phases dp
vstart=pagemtimes(QWP(pi/4),vstart); %QWP bias



dd=1/32;
phases(1,1,:)=-1+dd:dd:1; 
dphi=-1/4:1/4:1/2; %pi

outchannel=1;
figure(88);
clf
subplot(1,2,1);hold off;
sphases=squeeze(phases);
col=lines(length(dphi));
legends={};
for k=1:length(dphi)
    [vout, dangle]=EomSLM(vstart,phases*pi,dphi(k)*pi);
    int0=squeeze(abs(vout(outchannel,:,:)).^2/2);
    plot(sphases,int0,'Color',col(k,:))
    [~,ind]=min(int0);
    text(phases(ind),0.03,num2str(dphi(k)),"HorizontalAlignment","center","Color",'b')
    hold on
    legends{end+1}=['latpos = ' num2str(dphi(k)) ' pi'];
end
plot(sphases,0*sphases,'k')
text(min(phases),0.03,'lateral pos',"HorizontalAlignment","right","Color",'b')
ylabel('intensity')
xlabel('phase shift EOM / pi')
title('different positions')
legend(legends,'Location','north')

clear phases dphi
phases=-1/2:1/2:1; 
dphi(1,1,:)=-1/2:1/32:1/2; %pi

figure(88);subplot(1,2,2);hold off;
sphi=squeeze(dphi);
col=lines(length(dphi));
legends={};
for k=1:length(phases)
    [vout, dangle]=EomSLM(vstart,phases(k)*pi,dphi*pi);
    int0=squeeze(abs(vout(outchannel,:,:)).^2/2);
    plot(sphi,int0,'Color',col(k,:))
    [~,ind]=min(int0);
    text(dphi(ind),0.03,num2str(phases(k)),"HorizontalAlignment","center","Color",'b')        
    hold on
    legends{end+1}=['phase = ' num2str(phases(k)) ' pi'];
end
text(min(dphi),0.03,'EOM phase',"HorizontalAlignment","right","Color",'b')
plot(sphi,0*sphi,'k')
title('different EOM phases')
legend(legends,'Location','north')
xlabel('lateral position / pi')
ylabel('intensity')


 figure(89);hold off
 [vout, dangle,vplot]=EomSLM(vstart,pi/2,0);
 for k=1:length(vplot)
     plotpolstate(vplot(k).state,2*k,0)
     text(2*k-1,2,vplot(k).label)
     hold on
 end
axis equal
xlim([0 16])


function [vout, dangle,vplot]=EomSLM(vstart,phases,dphi) %EOM at 45°
if nargin<3
    dphi=0;
end
vplot(1).state=vstart; vplot(1).label='in';
vphase=pagemtimes(waveplate(phases,pi/4),vstart); %EOM 45  
vplot(2).state=vphase; vplot(2).label='EOM 45°';
vphase=pagemtimes(HWP(pi/8),vphase); %HWP
vplot(3).state=vphase; vplot(3).label='HWP 22.5°';
vrot=pagemtimes(addphase(-dphi),pagemtimes(HWP45,vphase)); %SLM+lat pos
vplot(4).state=vrot; vplot(4).label='SLM+latpos';
vphase=pagemtimes(addphase(dphi),vphase); %lat pos
vplot(5).state=vphase; vplot(5).label='only latpos';
vsum=(vrot+vphase);
vplot(6).state=vsum; vplot(6).label='sum';
vout=pagemtimes(polarizer(0),vsum); %polarizer
vplot(7).state=vout; vplot(7).label='polarizer';
outchannel=1;
dangle=mod(angle(vrot(outchannel,:))-angle(vphase(outchannel,:))+pi,2*pi)-pi;
end

function out=waveplate(phaseretardation,rotationangle)
out=exp(-1i*phaseretardation/2).*...
    [cos(rotationangle).^2+exp(1i*phaseretardation)*sin(rotationangle).^2 (1-exp(1i*phaseretardation))*cos(rotationangle).*sin(rotationangle)
    (1-exp(1i*phaseretardation))* cos(rotationangle).*sin(rotationangle) sin(rotationangle).^2+exp(1i*phaseretardation)*cos(rotationangle).^2]; %checked, correct
end

function out=QWP(rotationangle)
out=waveplate(pi/2,rotationangle);
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
out=waveplate(phase,0); 
% out=waveplate(2*phase,0);%Why factor of two?
end
function out=addphase(phase)
out=exp(1i*phase);
end

function out=linear(angle)
v=[1;0];
out=1i*HWP(angle/2)*v; % I changed it to 2 so it matches, not tested.
end

function out=polarizer(angle)
out=[cos(angle)^2 cos(angle)*sin(angle)
   cos(angle)*sin(angle) sin(angle)^2];
end


function plotpolstate(vin,dx,dy)
dd=0.1;
phases=0:dd:2*pi+dd;
p1=angle(vin(1));
p2=angle(vin(2));
i1=abs(vin(1));
i2=abs(vin(2));
plot(i1*sin(phases+p1)+dx,i2*sin(phases+p2)+dy); hold on
% 
plot(i1*sin(phases(1)+p1)+dx,i2*sin(phases(1)+p2)+dy,'b*'); hold on
plot(i1*sin(phases(2)+p1)+dx,i2*sin(phases(2)+p2)+dy,'bo'); hold on

end