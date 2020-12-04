
pixs=100;
dz=10;

%% BEADS: calibrate sx
figure(88);

sx=g.locData.loc.PSFxnm/pixs;
sy=g.locData.loc.PSFynm/pixs;
frame=g.locData.loc.frame;
z=frame*dz;
plot(frame,sx,'.',frame,sy,'.')

ds=sx.^2-sy.^2;
figure(89)
plot(frame,ds,'.')

range=[45,95];

inr=frame>range(1)&frame<range(2);

figure(89)
hold off
plot(ds(inr),z(inr),'.')

calbead=fit(ds(inr),z(inr),'poly4');
hold on
plot(ds(inr),calbead(ds(inr)))

maxrangeds=[min(ds(inr)) max(ds(inr))];

%% data: calculate z
outz=1000;
sx=g.locData.loc.PSFxnm/pixs;
sy=g.locData.loc.PSFynm/pixs;
ds=sx.^2-sy.^2;
z=calbead(ds);
g.locData.loc.znm=z-750;

outofrange=ds<maxrangeds(1) | ds>maxrange(2);
g.locData.loc.znm(outofrange)=outz;

%% data: average x,y
trafo=g.locData.files.file.transformation;
inref=trafo.getRef(g.locData.loc.xnm,g.locData.loc.ynm);
cr=[g.locData.loc.xnm(inref),g.locData.loc.ynm(inref)];
ct=[g.locData.loc.xnm(~inref),g.locData.loc.ynm(~inref)];
ctt=trafo.transformToReference(2,ct/pixs)*pixs;
lr.x=g.locData.loc.xnm;lr.y=g.locData.loc.ynm;lr.frame=g.locData.loc.frame;
ft=g.locData.loc.frame(~inref);
locpt=g.locData.loc.locprecnm(~inref);
lt.x=ctt(:,1);lt.y=ctt(:,2);lt.frame=ft;
[iA,iB]=matchlocsall(lr,lt,0,0,500);

normw=1./g.locData.loc.locprecnm(iA).^2+1./locpt(iB).^2;
xav=(lr.x(iA)./g.locData.loc.locprecnm(iA).^2+lt.x(iB)./locpt(iB).^2)./normw;
yav=(lr.y(iA)./g.locData.loc.locprecnm(iA).^2+lt.y(iB)./locpt(iB).^2)./normw;

figure(91)
plot(lr.y(iA),lt.y(iB),'.')


% plot(g.locData.loc.xnm(inref),g.locData.loc.ynm(inref),'.',g.locData.loc.xnm(~inref),g.locData.loc.ynm(~inref),'.')
figure(92)
plot(g.locData.loc.xnm(inref),g.locData.loc.ynm(inref),'.',cn(:,1),cn(:,2),'.',xav,yav,'.')


ind=true(size(g.locData.loc.xnm));
ind(iA)=false;
g.locData.removelocs(ind)
g.locData.loc.xnm=xav; 
g.locData.loc.ynm=yav;

g.locData.regroup;
