%Jones martrices
v=[0;1];
vlin45=linear(-pi/4);  
% vlin45=linear(-pi/4-0.25); 
eomqwp=1;

vstart=vlin45;
vstart=v;
clear phases dp
% vstart=pagemtimes(QWP(pi/4),vstart); %QWP bias

% figure(22);
% phases(1,1,:)=0:0.1:pi;
% hold off
% out=pagemtimes(waveplate(pi*0.7,pi/4),v)+1*v;
% plotpolstate(out,0,0)
% axis equal
% asfadf


dd=1/32;
phases(1,1,:)=-1+dd:dd:1; 
dphi=-1/4:1/4:1/2; %pi
% dphi=[-0.3 0 0.3];

outchannel=1;
figure(88);
clf
subplot(1,2,1);hold off;
sphases=squeeze(phases);
col=lines(length(dphi));
legends={};
for k=1:length(dphi)
    [vout, dangle]=EomSLM(vstart,phases*pi,dphi(k)*pi,eomqwp);
    int0=squeeze(abs(vout(outchannel,:,:)).^2/2);
    plot(sphases,int0,'Color',col(k,:))
    [~,ind]=min(int0);
    text(phases(ind),0.03,num2str(dphi(k)),"HorizontalAlignment","center","Color",'b')
    hold on
    legends{end+1}=['latpos = ' num2str(dphi(k)) ' pi, C = ' num2str(min(int0)/max(int0)*100,2) '%'];
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
    [vout, dangle]=EomSLM(vstart,phases(k)*pi,dphi*pi,eomqwp);
    int0=squeeze(abs(vout(outchannel,:,:)).^2/2);
    plot(sphi,int0,'Color',col(k,:))
    [~,ind]=min(int0);
    text(dphi(ind),0.03,num2str(phases(k)),"HorizontalAlignment","center","Color",'b')        
    hold on
    legends{end+1}=['phase = ' num2str(phases(k)) ' pi, C = ' num2str(min(int0)/max(int0),2)];
end
text(min(dphi),0.03,'EOM phase',"HorizontalAlignment","right","Color",'b')
plot(sphi,0*sphi,'k')
title('different EOM phases')
legend(legends,'Location','north')
xlabel('lateral position / pi')
ylabel('intensity')


 figure(85);hold off
 latpos=0;
 eomphase=0.15*pi;
 useqwp=0;
 qwpEOM=1;

 dx=1.7;
 % vstart=[0;1];
 vstart=vlin45;
 vstart=v;
 txtqwp='No QWP bias. ';
 if useqwp
    vstart=pagemtimes(QWP(pi/4),vstart); %QWP bias
     txtqwp='QWP bias. ';
 end

 [vout, dangle,vplot]=EomSLM(vstart,eomphase,0,qwpEOM);

 for k=1:length(vplot)
     plotpolstate(vplot(k).state,dx*k,0)
     text(dx*k-dx/2,dx,vplot(k).label,"HorizontalAlignment","center")
     hold on
 end
axis equal
xlim([0 dx*k+dx])
title([txtqwp 'phase ' num2str(eomphase/pi) '\pi, lat pos ' num2str(latpos/pi) '\pi'])




function [vout, dangle,vplot]=EomSLM(vstart,phases,dphi,useqwp) %EOM at 45°
if nargin <4
    useqwp=false;
end
if nargin<3
    dphi=0;
end
vplot(1).state=vstart; vplot(1).label='in 45°';
vphase=vstart;

if useqwp
    % vphase=pagemtimes(QWP(0),vphase); 
    vphase=pagemtimes(QWP(pi/4),vphase); 
    % vphase=pagemtimes(waveplate(pi/4,0),vphase); %EOM 0  
    % vphase=pagemtimes(waveplate(pi/4,0),vphase); 
    vplot(2).state=vphase; vplot(2).label='EOM QWP 0°';
else
    vplot(2).state=vphase; vplot(2).label='-';
end

vphase=pagemtimes(waveplate(phases,0),vphase); %EOM 0  
vplot(3).state=vphase; vplot(3).label='EOM 0°';


vphase=pagemtimes(HWP(-0.15),vphase); %HWP to compensate for non-ideal slm
% vphase=pagemtimes(QWP(-0.2),vphase); %HWP

% vplot(4).state=vphase; vplot(4).label='HWP 22.5°';
vplot(4).state=vphase; vplot(4).label='-';

vrot=pagemtimes(addphase(-dphi),pagemtimes(waveplate(pi*555/640,33.5/180*pi),vphase)); %realistic SLM: angle 33.5, retardance= pi*555/640
% vrot=pagemtimes(addphase(-dphi),pagemtimes(HWP45,vphase)); %SLM+lat pos
% vrot=pagemtimes(addphase(-dphi),pagemtimes(waveplate(pi/2-0*pi/8,pi/4),vphase)); %SLM broken+lat pos
vplot(5).state=vrot; vplot(5).label='SLM+latpos';
vphase=pagemtimes(addphase(dphi),vphase); %lat pos
% vphase=pagemtimes(waveplate(pi*555/640*0,0),vphase); % realistic SLM: off pixel: optical axis vertical
vplot(6).state=vphase; vplot(6).label='only latpos';
vsum=(vrot+vphase);
vplot(7).state=vsum; vplot(7).label='sum';
vout=pagemtimes(polarizer(0),vsum); %polarizer
vplot(8).state=vout; vplot(8).label='polarizer';
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
if nargin<3
    dy=0;
end
if nargin <2
    dx=0;
end
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