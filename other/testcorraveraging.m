% create ring

dp=pi/180/2; %.5 pixel size in rad for polar angle 
sizegaussdeg= 10;  % x degree
sizegauss1=sizegaussdeg/180*pi/dp;
sizegauss2=sizegauss1;
% dp=2*dp
np=0:dp:pi/4;
npdegh=np(1:end-1)*180/pi;
numparticles=1000;
maxiter=10;
plabel=0.5;
%% make NPCs
if 0
    phierr=pi/32;
    for k=numparticles:-1:1
        phis{k}=simulateNPC(phierr,plabel);
    end

%get from CSV
else
    file='/Users/jonasries/Library/Containers/nz.co.pixeleyes.AutoMounter/Data/Mounts/ries.embl.de/SMB/ries/users/Yu-Le/NUP96_hires/11_U2OS_GFPNanoBQAF647__lowerRing_sml.csv';
    tc=readtable(file);
    classes=unique(tc.class);
    [~,indsort]=sort(rand(length(classes),1));
    classes=classes(indsort);
    numparticles=length(classes);
   
    for k=numparticles:-1:1
        ind=tc.class==classes(k);
        [phis{k},radius{k}]=cart2pol(tc.xnm(ind)-500,tc.ynm(ind)-500);
        phis{k}=phis{k};%+rand*pi/4;
        numlocsparticle(k)=sum(ind);
    end
end

%
[classsort,indsort]=sort(numlocsparticle,'descend');
phis=phis(indsort);
% minlocs=20;
% indgood=numlocsparticle>=minlocs;
% phis=phis(indgood);
%     numparticles=length(phis);
%% align
phiout=alignhalf(phis,np,sizegauss1);
figure(299);href=histcounts(mod(phiout,pi/4), np);plot(npdegh,href);title('pairwise alignment');

%% align iteratively
href=sharpenref(href);

phisrot=phis;
htar=0*href;
for iter=1:maxiter
    for k=1:numparticles
         hphit=histcounts(mod(phisrot{k},pi/4), np);
         drot=getrotangle(href,hphit,np,sizegauss2);
         phisrot{k}=phisrot{k}-drot;
         hphitr=histcounts(mod(phisrot{k},pi/4), np);
         htar=htar+hphitr;
    end
    figure(298);plot(npdegh,htar);title('iterative refinement');drawnow
    href=htar;
    href=sharpenref(href);
end

%% simple cross-correalation
% [phir,rotangle]=simulateNPC(phierr);
hphi=histcounts(mod(phis{1},pi/4), np);
for k=2:numparticles
%     [phit,rotanglet]=simulateNPC(phierr);
%     drotangle_gt=mod(rotanglet-rotangle,pi/4);
    hphit=histcounts(mod(phis{k},pi/4), np);
    hphiref=hphi/k;
%     hphiref(hphiref<0.25)=0;
%     hphiref(hphiref>2)=2;
%     hphiref=(hphi/k)>0.5;
%     hphiref=sqrt(hphi)/k;
    drot=getrotangle(hphiref,hphit,np,sizegauss1);
    hphitr=histcounts(mod(phis{k}-drot,pi/4), np);
    hphi=hphi+hphitr;
%     getrotangle(hphi,hphitr,np)
end

figure(188);plot(npdegh,hphi/k);title('simple cross-correlation');

function phiout=alignhalf(phis,np,sizegauss)
s=length(phis);
if s(1)==1
    phiout=phis{1};
    return
elseif s(1)==2
    phi1=phis{1};
    phi2=phis{2};
else
    mp=floor(s(1)/2);
    phi1=alignhalf(phis(1:mp),np,sizegauss);
    phi2=alignhalf(phis(mp+1:end),np,sizegauss);   
end
    hr=histcounts(mod(phi1,pi/4), np);
    ht=histcounts(mod(phi2,pi/4), np);
    drot=getrotangle(hr,ht,np,sizegauss);
    phiout=[phi1; phi2-drot];
end

function [phio,rotangle]=simulateNPC(phierr,plabel)
dpair=0.25;
rotangle=rand*(pi/4);
phi=0:pi/4:7/4*pi;
phi=phi'+rotangle;
phi16=[phi ;phi+dpair];
phi16(rand(length(phi16),1)>plabel)=[];
phio=phi16+randn(length(phi16),1)*phierr;
end

function drot=getrotangle(href,htar,np,sizegauss)

href=sqrt(href);
htar=sqrt(htar);
nump=length(href);
numphalf=round(nump/2);
 cr=xcorrangle(href,htar);
    crf=smoothdata([cr,cr],'gaussian',sizegauss);
    crfs=[crf(nump+1:nump+numphalf) crf(numphalf+1:nump)];
    [~,indmax]=max(crfs);
    drot=np(indmax);
end

function hout=sharpenref(hin)
% hout=hin.^2;
hout=hin;
% hout(hout<max(hout)/2)=0;
% hout=sqrt(hin);
% hout=hout-max(hout)/2;
% hout(hout<0)=0;
end