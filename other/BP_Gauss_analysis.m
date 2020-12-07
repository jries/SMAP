
% pixs=100;
dz=10;
zoff=750;

%% BEADS: calibrate sx
pixs=g.locData.files.file(1).info.cam_pixelsize_um(1)*1000;
loc=g.locData.getloc({'PSFxnm','PSFynm','frame'},'layer',1,'position','roi');
indg=loc.PSFxnm>0.1&loc.PSFynm>0.1;
sx=loc.PSFxnm(indg)/pixs;
sy=loc.PSFynm(indg)/pixs;


frame=loc.frame(indg);
z=frame*dz-zoff;
figure(88);
plot(frame,sx,'+',frame,sy,'.')

ds=sx.^2-sy.^2;
figure(89)
plot(frame,ds,'.')

range=[44,100];

inr=frame>range(1)&frame<range(2);

figure(89)
hold off
plot(ds(inr),z(inr),'.')

calbead=fit(ds(inr),z(inr),'smoothingspline','SmoothingParam',.1);
hold on
nds=min(ds(inr)):0.1:max(ds(inr));
plot(nds,calbead(nds))

maxrangeds=[min(ds(inr)) max(ds(inr))];
%%
pfad=[fileparts(g.locData.files.file(1).name) filesep 'beadcal_Gauss_smoothingspline'];
save(pfad,'calbead','maxrangeds')

%% data: calculate z
pixs=g.locData.files.file(1).info.cam_pixelsize_um(1)*1000;
outz=1000;
sx=g.locData.loc.PSFxnm/pixs;
sy=g.locData.loc.PSFynm/pixs;
ds=sx.^2-sy.^2;
z=calbead(ds);
g.locData.loc.znm=z;

outofrange=ds<maxrangeds(1) | ds>maxrangeds(2) | sx==0 | sy==0;
g.locData.loc.znm(outofrange)=outz;
g.locData.regroup;


%% data: average x,y
%use combine channels plugin 