function [img,beads]=validatemodel4Pi(PSF,ph,titlet)
if nargin<3
    titlet='results';
end
if isfield(ph,'Nfree')
    Nfree=ph.Nfree;
else
    Nfree=false;
end
if isfield(ph,'xyfree')
    xyfree=ph.xyfree;
else
    xyfree=false;
end
ph.isglobalfit=true;

[beads,ph]=images2beads_globalfitN(ph); 

[imstack,fn,dxy]=bead2stack(beads);
img.imstack=imstack;
sim=size(imstack);
imsqueeze=reshape(imstack,sim(1),sim(2),[],sim(end));

if isfield(PSF,'normf')
for k=2:size(imsqueeze,4)
    imsqueeze(:,:,:,k)=imsqueeze(:,:,:,k)/PSF.normf(k);
end
end

dTAll=reshape(dxy,size(dxy,1),sim(end),[]);
img.dTAll=dTAll;
if Nfree
    shared=[1,1,0,1,1,1];
    indN=3:6;
    indz=8;
    indp=9;
else
    shared=[1,1,1,1,1,1];
    indN=3;
    indz=5;
    indp=6;
end
if xyfree
    shared(1:2)=0;
    indN=indN+6;
    indz=indz+6;
    indp=indp+6;
end
imstacksq=imsqueeze(ph.rangeh, ph.rangeh, :, :);

z0a=ph.zstart;

%prefit with free N
% sharedh=shared;

% iterationsh=100;
% sharedh(3)=0;
% tic
% z0=[-5 ];
% 
% P0 = mleFit_LM_4Pi(single(imstacksq(:, :, :, :)),uint32(ones(size(shared))),iterationsh,single(PSF.Ispline), single(PSF.Aspline),single(PSF.Bspline),single(dTAll),single(ph.phi0),0);
% 
% z0a=P0(:,end-2)-(size(PSF.I,3)/2);
% p0=mod(P0(:,end-1),2*pi);

% Pi = mleFit_LM_4Pi(single(imstacksq(3:end-2, 3:end-2, :, :)),uint32(sharedh),iterationsh,single(PSF.Ispline), single(PSF.Aspline),single(PSF.Bspline),single(dTAll),single(ph.phi0),z0a0,p00);
% toc
% tic
% Pi2 = mleFit_LM_4Pi(single(imstacksq(:, :, :, :)),uint32(sharedh),iterationsh,single(PSF.Ispline), single(PSF.Aspline),single(PSF.Bspline),single(dTAll),single(ph.phi0),z0);
% toc

% z0a=Pi(:,end-2)-(size(PSF.I,3)/2);
% p0=Pi(:,end-1);
% Nnotlinked=reshape(Pi(:,indN(1):indN(1)+3),[],sim(4),4);

% z0a=0;p0=0;
% z0a=z0;
p0a=[0 pi]+pi/4;
iterations=100;
% sharedh(3)=1;
[P,CRLB1 LL] = mleFit_LM_4Pi(single(imstacksq(:, :, :, :)),uint32(shared),iterations,single(PSF.Ispline), single(PSF.Aspline),single(PSF.Bspline),single(dTAll),single(ph.phi0),z0a,p0a);

% P=P0;
% [P,CRLB1 LL] = CPUmleFit_LM_MultiChannel_4pi(single(imstacksq(:, :, :, :)),uint32(shared),iterations,single(PSF.Ispline), single(PSF.Aspline),single(PSF.Bspline),single(dTAll),single(ph.phi0),z0);

img.imstacksq=imstacksq;
img.sim=sim;
% img.fit.Nnotlinked=Nnotlinked;
img.fit.P=P;
img.fit.CRLB=CRLB1;
img.fit.PSF=PSF;
img.fit.shared=shared;
%now unlink x, y to see if there is shift
% shared(1:2)=0;
% [Pu,CRLB1 LL] = CPUmleFit_LM_MultiChannel_4pi(single(imstacksq(:, :, :, :)),uint32(shared),iterations,single(PSF.Ispline), single(PSF.Aspline),single(PSF.Bspline),single(dTAll),single(phi0),z0);
% dx21=Pu(:,2)-Pu(:,1);
 
%collect fitted parameters
phase=mod(reshape(P(:,indp),[],sim(4)),2*pi);
zphase=phase/2/PSF.frequency*ph.dz;
zastig=reshape(P(:,indz),[],sim(4))*ph.dz;
xfit=reshape(P(:,1),[],sim(4));
yfit=reshape(P(:,2),[],sim(4));
N=reshape(P(:,indN),[],sim(4),length(indN));


LLi=reshape(LL,[],sim(4));
iter=reshape(P(:,end),[],sim(4));

%XXXX find z0!
z_phi = reshape(z_from_phi_JR(P(:, indz), phase(:), PSF.frequency, ceil(sim(3)/2)-.7),[],sim(4))*ph.dz;

%plot results of validation
tab=(uitab(ph.tabgroup,'Title',['r_' titlet]));
tgr=uitabgroup(tab);
ax=axes(uitab(tgr,'Title','z_astig'));
plot(ax,zastig)
xlabel(ax,'frame')
ylabel(ax,'z_astig')
ax=axes(uitab(tgr,'Title','phase'));
plot(ax,phase)
xlabel(ax,'frame')
ylabel(ax,'phase')
ax=axes(uitab(tgr,'Title','phase(z_a)'));
plot(ax,zastig,zphase)
xlabel(ax,'z_astig')
ylabel(ax,'z_phase')

ax=axes(uitab(tgr,'Title','z_phase'));
plot(ax,z_phi)
xlabel(ax,'frame')
ylabel(ax,'z_phi')
ax=axes(uitab(tgr,'Title','dz_phase'));
frame=1:length(z_phi);
dzp=z_phi-frame'*ph.dz;
mz=ceil(size(z_phi,1)/2);
dzp=dzp-dzp(mz,:);
plot(ax,dzp)
xlabel(ax,'frame')
ylabel(ax,'z_phi-dz*frame (nm)')

ax=axes(uitab(tgr,'Title','x,y'));
plot(ax,xfit,yfit,'+')
xlabel(ax,'x')
ylabel(ax,'y')
ax=axes(uitab(tgr,'Title','x(z)'));
hold(ax,'off')
plot(ax,zastig,xfit)
hold(ax, 'on')
xlabel(ax,'z_astig')
ylabel(ax,'x')

ax=axes(uitab(tgr,'Title','LL'));

histogram(ax,LL/sim(1)^2)
title(ax,median(LL)/sim(1)^2)

if Nfree
    ax=axes(uitab(tgr,'Title','N'));
%     if Nfree
        Np=squeeze(mean(N(28:34,:,:),1));
%     else
%         Np=squeeze(mean(Nnotlinked(28:34,:,:),1));
% end 
        
    plot(ax,Np./Np(:,4))
    xlabel(ax,'bead')
    ylabel(ax,'Nfit')
end

% calculate residuals for each bead. Should it beased on average x,y,z in center?

end