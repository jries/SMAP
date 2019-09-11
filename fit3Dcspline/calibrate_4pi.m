function calibrate_4pi(p)
%add directory with fitter to path
fit4pidir=strrep(pwd,'SMAP',['ries-private' filesep 'PSF4Pi']);
if exist(fit4pidir,'file')
    addpath(fit4pidir);
end
ph=p;
f=figure(144);
tg=uitabgroup(f);
t1=uitab(tg,'Title','prefit');
tgprefit=uitabgroup(t1);
ph.tabgroup= tgprefit;
ph.isglobalfit=false;
ph.outputfile={};
% ph.normalize=false;

% get parameters for cutting out and mirroring raw image files
if ~isempty(p.settingsfile4pi) && exist(p.settingsfile4pi,'file')
    settings_3D=readstruct(p.settingsfile4pi); %later to settings, specify path in gui  
    settings_3D.file=p.settingsfile4pi;
    ph.settings_3D=settings_3D;
end

ph.zstart= [-1 0 1]*30;
% ph.zstart= 0;
%segment beads
 p.status.String=['load files'];drawnow
[beads,ph]=images2beads_globalfit(ph); %later extend for transformN, dont use two versions of images2beads 
 p.status.String=['register beads in x,y,z'];drawnow
%register beads channel wise
filenumbers=[beads(:).filenumber];
infile=filenumbers==1;
[allPSFs,shiftedstack,corrout]=PSFcorrelation(beads(infile),ph);

%get frequency and phases
tab=(uitab(tgprefit,'Title','frequency'));ph.ax=axes(tab);
 p.status.String=['Get phases'];drawnow
[phaseh,ph.frequency]=getphaseshifts(allPSFs,ph.ax,ph);
phaseshifts=[phaseh(1) phaseh(2) phaseh(1)+pi phaseh(2)+pi]; 
phaseshifts=phaseshifts-phaseshifts(1)-pi;


%first step: align channels via CC, calculate transform, find corresponding
%beads. As before.

%align four quadrants using IABfrom4PIPSFfit
%as averages are calculated channel-wise, we have to make sure only
%corresponding beads / regions are taken into account. Otherwise intensities totally off. This is now taken care of in PSFcorelation  

%alignment in z: not good. individual beads are already shifted, but all
%channels need to be shifted the same way. Now that same beads under
%consideration, look at mean z shift. This works, made it much mor robust.
%also, consider doing some normalization, so not intensity but bead counts
%in the average
 p.status.String=['make spline model'];drawnow
PSF=IABfrom4PiPSFfit(allPSFs, phaseshifts(2),ph.frequency,9,50,corrout.zshift0);
[PSFspl,globalnorm]=makeIABspline(PSF.I,PSF.A,PSF.B,p);
PSF=copyfields(PSF,PSFspl);

for k=1:corrout.numchannels
    corrout.beadtrue{k}(:,1)= corrout.beadtrue{k}(:,1)-PSF.dx(k);%shift2(k,2); %out.dx
    corrout.beadtrue{k}(:,2)=corrout.beadtrue{k}(:,2)-PSF.dy(k);%shift2(k,1); %out.dy
    corrout.beadtrue{k}(:,3)=corrout.beadtrue{k}(:,3)-PSF.dz(k);%shift2(k,3)+corrout.zshift0(k); %  out.dz %not sure about sign %needed?
end
 p.status.String=['Calculate transformation'];drawnow
ph.transformation=make4PiTransform(corrout.beadtrue,ph);
out.transformation=ph.transformation;

%now: validation and plotting of graphs
%do fitting for testing
fitroi=13;
sim=size(allPSFs);
mp=floor((sim(1)-1)/2)+1;
mpz=floor((sim(3)-1)/2)+1;
droi=floor((fitroi-1)/2);
ph.rangeh=mp-droi:mp+droi;
ph.phi0=phaseshifts;



%plot PSF
plotI(:,:,:,1)=PSF.I/globalnorm;plotI(:,:,:,2)=PSF.A/globalnorm;plotI(:,:,:,3)=PSF.B/globalnorm;plotI(:,:,:,4)=PSF.PSF(:,:,:,1)/globalnorm;
plotI(:,:,:,5)=allPSFs(:,:,:,1)/globalnorm;
tab=(uitab(tgprefit,'Title','IAB'));imageslicer(plotI,'Parent',tab)

for k=1:4
plotR(:,:,:,k)=(allPSFs(:,:,:,k)-PSF.PSF(:,:,:,k))/globalnorm;
end
plotR(:,:,1,:)=[];plotR(:,:,end,:)=[];
tab=(uitab(tgprefit,'Title','residuals'));
tgres=uitabgroup(tab);
tabres=(uitab(tgres,'Title','residuals all'));
imageslicer(plotR,'Parent',tabres)
    
%residuals for all beads in shiftedstack
%try for bead 1
rxy=ph.rangeh;rz=mpz-floor(p.zcorrframes/2):mpz+floor(p.zcorrframes/2);

for bead=1:size(shiftedstack{1},4)
    for ch=1:4
        ratiod(:,:,:,ch)=shiftedstack{ch}(rxy,rxy,rz,bead)./PSF.PSF(rxy,rxy,rz,ch);
        int1=shiftedstack{ch}(rxy,rxy,rz,bead);int2=PSF.PSF(rxy,rxy,rz,ch);
         ratio(ch)=int1(:)\int2(:);
    end
    normhd=median(ratiod(:),'omitnan');
    normh=1/mean(ratio);
    resh=[];
    for ch=1:4
        resh=vertcat(resh,shiftedstack{ch}(:,:,:,bead)-normh*PSF.PSF(:,:,:,ch));
        ressmallh=shiftedstack{ch}(rxy,rxy,rz,bead)-normh*PSF.PSF(rxy,rxy,rz,ch);
        ressum(bead,ch)=sqrt(sum(ressmallh(:).^2,'omitnan'))/globalnorm;
    end
    residuals(:,:,:,bead)=resh;
end
tab=(uitab(tgres,'Title','res beads'));imageslicer(residuals,'Parent',tab)
tab=(uitab(tgres,'Title','ressum'));ax=axes(tab);plot(ax,ressum);

 p.status.String=['Validate by fitting'];drawnow
valfit=validatemodel(PSF,ph,'fit');
%fit calibrations stack
% shared=[0,0,1,1,1,1];
% z0=ph.zstart;
% dTAll=zeros(6,4,size(allPSFs,3),'single');
% iterations=50;
% % imstack=allPSFs(ph.rangeh, ph.rangeh, :, :)*10000;
% imstack=PSF.PSF(ph.rangeh, ph.rangeh, :, :)*10000;
% [Pc,CRLB1 LL] = mleFit_LM_4Pi(single(imstack(:, :, :, :)),uint32(shared),int32(iterations),single(PSF.Ispline), single(PSF.Aspline),single(PSF.Bspline),single(dTAll),single(ph.phi0),single(z0));
% nanmean(Pc(:,1:8),1)-droi+1 % if this is not all the same -> PSFs in channels not well aligned. 

%cut out corresponding beads based on transform to mimick normal fitting



% % old approach: align PSFs. For comparison
% sstack=size(beads(1).stack.image);
% fw=round((ph.zcorrframes-1)/2);
% framerange=round(max(1,sstack(3)/2-2*fw):min(sstack(3)/2+2*fw,sstack(3)));
% [~,PSFaligned,shift2,indgood]=registerPSF3D_g(allPSFs,[],...
%     struct('framerange',framerange,'removeoutliers',false,'alignz',false,...
%     'zshiftf0',+corrout.zshift0),{},corrout.filenumber);

% %make IAB model from average PSFs
% 
% [Io,Ao,Bo,PSFo]=make4Pimodel(PSFaligned,phaseshifts,ph.frequency,p);
% % PSFo.normf=PSF.normf;
% img=validatemodel(PSFo,ph,'corr');
% 
%     plotI(:,:,:,1)=Io;plotI(:,:,:,2)=Ao;plotI(:,:,:,3)=Bo;plotI(:,:,:,4)=PSFaligned(:,:,:,1);
%     tab=(uitab(tgprefit,'Title','IAB'));imageslicer(plotI,'Parent',tab)
    
    
    %plot
    tab=(uitab(tgprefit,'Title','PSFaligned'));
    imageslicer(cat(4,shiftedstack{1},shiftedstack{2}),'Parent',tab)
    tab=(uitab(tgprefit,'Title','Zprofile'));
    simh=size(shiftedstack{1});
    mp=ceil((simh(1)+1)/2);
    ax=axes(tab);
    hold (ax,'off')
    for k=1:4
        dz=1*(k-1);
        profilez=squeeze(shiftedstack{k}(mp,mp,:,:));
        profilezm=squeeze(allPSFs(mp,mp,:,k));
        profilezf=squeeze(PSF.PSF(mp,mp,:,k));
        plot(ax,profilez./max(profilez)+dz)
        hold (ax,'on')
        plot(ax,profilezm/max(profilezm)+dz,'k')
        plot(ax,profilezf/max(profilezf)+dz,'rx')
    end
    norm=max(profilezm);
    for k=1:2
       
        profilez=(squeeze(mean(shiftedstack{k}(mp,mp,:,:),4))+squeeze(mean(shiftedstack{k+2}(mp,mp,:,:),4)));
        plot(ax,profilez/norm*4)
    end

out.cal4pi.coeff=PSF;
out.cal4pi.dz=ph.dz;
out.cal4pi.x0=ceil((ph.ROIxy+1)/2);
out.cal4pi.z0=ceil((size(PSF.Aspline,3)+2)/2);
out.cal4pi.transformation=out.transformation;
out.cal4pi.settings3D=ph.settings_3D;
out.Xrange=ph.xrange;
out.Yrange=ph.yrange;
out.EMon=ph.emgain;
parameters=rmfield(ph,{'tabgroup','ax','status','fileax','smappos'});
out.parameters=parameters;
  
p.status.String='save calibration';drawnow
if ~isempty(p.outputfile)
        save(p.outputfile,'-struct','out');
end

 p.status.String=['Calibration done.'];drawnow
 
 
 

% 
% % OOOOOOOO now testing global fit of all beads

% pass on start parameters (e.g. previous fit parameters from validation)
% optimize also for internal consistency
% try with global norm
% 
% ph.isglobalfit=true;
% [beads,ph]=images2beads_globalfitN(ph); %get global bead stacks
% [imstack,fn,dxy]=bead2stack(beads);
% PSF.globalnorm=globalnorm;

startp=averagefit4Pi(valfit);
out=IABfrom4PiPSFfitmany(valfit,startp, phaseshifts(2),ph.frequency,9,30);

% OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO 
% now testing iterative fitting and alignemnt of beads with
% averaging IABs
ph.isglobalfit=true;
[beads,ph]=images2beads_globalfitN(ph); %get global bead stacks
[imstack,fn,dxy]=bead2stack(beads);

valfit.imstack=imstack;
sim=size(imstack);
imsqueeze=reshape(imstack,sim(1),sim(2),[],sim(end));
iterations=5
% it seems that with every iteration it gets worse (x,y registration etc)
for iter=1:iterations
% for k=2:size(imsqueeze,4)
%     imsqueeze(:,:,:,k)=imsqueeze(:,:,:,k)/PSF.normf(k);
% end

dTAll=reshape(dxy,size(dxy,1),sim(end),[]);
valfit.dTAll=dTAll;
shared=[0,0,0,1,1,1];
imstacksq=imsqueeze(ph.rangeh, ph.rangeh, :, :);
iterations=50;
z0=ph.zstart;
[P,CRLB1 LL] = mleFit_LM_4Pi(single(imstacksq(:, :, :, :)),uint32(shared),iterations,single(PSF.Ispline), single(PSF.Aspline),single(PSF.Bspline),single(dTAll),single(ph.phi0),z0);

% [P,CRLB1 LL] = CPUmleFit_LM_MultiChannel_4pi(single(imstacksq(:, :, :, :)),uint32(shared),iterations,single(PSF.Ispline), single(PSF.Aspline),single(PSF.Bspline),single(dTAll),single(ph.phi0),z0);

 P=double(P);CRLB1=double(CRLB1);
%collect fitted parameters

for k=1:size(CRLB1,2)
    Pr(:,k,:)=reshape(P(:,k),[],valfit.sim(4));
    Cr(:,k,:)=reshape(CRLB1(:,k),[],valfit.sim(4));
end
Pr(:,k+1,:)=reshape(P(:,k+1),[],valfit.sim(4)); %iterations, not in crlb


xfit=Pr(:,1:4,:);dx=Cr(:,1:4,:);
yfit=Pr(:,5:8,:);dy=Cr(:,5:8,:);
Nfit=(Pr(:,9:12,:));
Bg=(Pr(:,13,:));
phase=mod(Pr(:,15,:),2*pi);
zastigf=squeeze(Pr(:,14,:));
%determine average positions in small window
zwindow=5; %+/- zwin
zrange=mpz-zwindow:mpz+zwindow;
numbeads=sim(4);


xn=1:size(valfit.imstack,1);yn=1:size(valfit.imstack,2);zn=1:size(valfit.imstack,3);
[Xq,Yq,Zq]=meshgrid(yn,xn,zn);

imstackaligned=imstack*0;
for k=numbeads:-1:1
    phasem(k)=cyclicaverage(phase(zrange,k),2*pi);
    zastigh=zastigf(zrange,k)-mpz;
    Xmat=horzcat(zastigh(:), zastigh(:)*0+1);
    linfit=Xmat\(zrange(:)-mpz);
    z0(k)=linfit(2);
%     x0(k,:)=squeeze(robustMean(xfit(zrange,k,:),1))-droi+1;
%     y0(k,:)=squeeze(robustMean(yfit(zrange,k,:),1))-droi+1;   
    x0(k,:)=squeeze(sum(xfit(zrange,:,k)./dx(zrange,:,k),1)./sum(1./dx(zrange,:,k),1))-droi+1;
    y0(k,:)=squeeze(sum(yfit(zrange,:,k)./dy(zrange,:,k),1)./sum(1./dy(zrange,:,k),1))-droi+1;
    
    Nhere(k,:)=squeeze(mean(Nfit(zrange,:,k),1));
    for c=1:size(valfit.imstack,5)
        imh=squeeze(imstack(:,:,:,k,c));
        xshift=-y0(k,c); %works empirically
        yshift=-x0(k,c);
        zshift=0; %shift IAB in z only
%         zshift=-z0(k);
        shiftedh=interp3(imh(:,:,:),Xq-xshift,Yq-yshift,Zq-zshift,'cubic',0);
        imstackaligned(:,:,:,k,c)=shiftedh;
    end
    [I,A,B]=make4Pimodel(squeeze(imstackaligned(:,:,:,k,:)),phaseshifts+phasem(k),ph.frequency,1./Nhere(k,:));
%        [I,A,B]=make4Pimodel(squeeze(imstackaligned(:,:,:,k,:)),phaseshifts+phasem(k),ph.frequency,PSF.normf);
    Is=interp3(I,Xq,Yq,Zq+z0(k),'cubic',0);
    As=interp3(A,Xq,Yq,Zq+z0(k),'cubic',0);
    Bs=interp3(B,Xq,Yq,Zq+z0(k),'cubic',0);
    
    Aa(:,:,:,k,:)=As;
    Ba(:,:,:,k,:)=Bs;
    Ia(:,:,:,k,:)=Is;
end
Am=squeeze(mean(Aa,4));
Bm=squeeze(mean(Ba,4));
Im=squeeze(mean(Ia,4));

% PSF=IABfrom4PiPSFfit(squeeze(sum(imstackaligned(:,:,:,:,:),4)), phaseshifts(2),ph.frequency,9,25,[0 0 0 0]);
[out,globalnorm2]=makeIABspline(Im,Am,Bm,ph);
PSFiter=copyfields(PSF,out);
valfit=validatemodel(PSFiter,ph,['fit' num2str(iter)]);
PSF=PSFiter;

end
% OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO



%first iteration
shared=[0 0 1 1 1 1]; %only link BG and photons to get true x,y
%phase shoudl be linked! comes from the cavity, the same in all quadrants
% different focus (=zposition) in four quadrants:destroys relationship
% between zastig and phase (phase stays constant, zastig changes). avoid!
% negelct here
%now fit with dT=0 to get directly the shift (avoid adding shifts)

dTAll0=valfit.dTAll*0;
[Pu,CRLB1 LL] = mleFit_LM_4Pi(single(valfit.imstacksq(:, :, :, :)),uint32(shared),iterations,single(PSF.Ispline), single(PSF.Aspline),single(PSF.Bspline),single(dTAll0),single(ph.phi0),z0);


for k=1:size(CRLB1,2)
    Pr(:,k,:)=reshape(Pu(:,k),[],valfit.sim(4));
    Cr(:,k,:)=reshape(CRLB1(:,k),[],valfit.sim(4));
end
Pr(:,k+1,:)=reshape(Pu(:,k+1),[],valfit.sim(4)); %iterations, not in crlb
% Pr=reshape(Pu,[],size(Pu,2),sim(4));
df=20;
frange=ceil((sim(3)-1)/2+1)+ (-df:df)';
mean(Pu(:,1:8),1)-droi+1


x=Pr(frange,1:4,:);dx=Cr(frange,1:4,:);
y=Pr(frange,5:8,:);dy=Cr(frange,5:8,:);
z=Pr(frange,11,:); dz=Cr(frange,11,:); dz(isnan(dz))=inf;
phase=mod(Pr(frange,12,:),2*pi); dphase=Cr(frange,12,:);

z_phase = (z_from_phi_JR((z), (phase), ph.frequency, ceil(sim(3)/2)-.7));

% xrm=squeeze(robustMean(x,1));
xwm=squeeze(sum(x./dx,1)./sum(1./dx,1))-droi+1;
ywm=squeeze(sum(y./dy,1)./sum(1./dy,1))-droi+1;

zr=z-frange;
zwm=squeeze(sum(zr./dz,1)./sum(1./dz,1));
[zwm,stdm,inlier]=(robustMean(zr,1));
zwm=squeeze(zwm);
stdm=squeeze(stdm);
indgood=abs(zwm)<200/ph.dz &stdm<2;
zpwm=squeeze(sum(z_phase./dphase,1)./sum(1./dphase,1));

% xm=squeeze(mean(x,1));
%determine crlb for weighing
%x,y, weighted average from central part
%z: z-frange, weighted average
%phase: cyclicaverage as in postfitting plugin. Maybe make robust? median?
%Or fit, get dhi from fit.


%phase vs zastig: very very well on line. Hardly any spread. Alignment in z
%sufficient? No need to adjust phase? 

%but this means that prefit already corrected for everything. But if phase
%not flat (imperfect alignment) across FoV then z would stay constant and
%phase would change. 

xn=1:size(valfit.imstack,1);yn=1:size(valfit.imstack,2);zn=1:size(valfit.imstack,3);
[Xq,Yq,Zq]=meshgrid(yn,xn,zn);

imstackaligned=valfit.imstack*0;
imstackalignedp=valfit.imstack*0;
%try also to align by phase
% AB, then average: zastig
%average then AB: zphase

% shift all images according to dT from fit
%imstackaligned: all beads, shifted by this, using average z_astig
%imstackalignedp: using z_phase
for k=1:size(valfit.imstack,4) %for all beads
    for c=1:size(valfit.imstack,5)
        imh=squeeze(valfit.imstack(:,:,:,k,c));
        xshift=-ywm(c,k); %works empirically
        yshift=-xwm(c,k);
        zshift=zwm(k);
        shiftedh=interp3(imh(:,:,:),Xq-xshift,Yq-yshift,Zq-zshift,'cubic',0);
        imstackaligned(:,:,:,k,c)=shiftedh;
        
        zshift=zpwm(k);
        shiftedh=interp3(imh(:,:,:),Xq-xshift,Yq-yshift,Zq-zshift,'cubic',0);
        imstackalignedp(:,:,:,k,c)=shiftedh;
    end
end
imstackaligned=imstackaligned(:,:,:,indgood,:);
imstackalignedp=imstackalignedp(:,:,:,indgood,:);

[imstackalignedn,factor]=normalizequadrants(imstackaligned);
%XXXX take into account factor when fitting!
%Put normalizequadrants into make4Pimodel? Or outside for global factor.
Iz=0;Bz=0;Az=0;
numbeads=size(imstackalignedn,4);
for k=1:numbeads %for each bead caluclate IAB, average those
    [Ih,Ah,Bh]=make4Pimodel(squeeze((imstackalignedn(:,:,:,k,:))),phaseshifts,ph.frequency);
    Iz=Ih/numbeads+Iz;Bz=Bh/numbeads+Bz;Az=Ah/numbeads+Az;
end
mp=ceil((size(Az,1)-1)/2);
dd=floor((ph.ROIxy-1)/2);
PSFz.Aspline=single(getsmoothspline(Az(mp-dd:mp+dd,mp-dd:mp+dd,:),ph));
PSFz.Bspline=single(getsmoothspline(Bz(mp-dd:mp+dd,mp-dd:mp+dd,:),ph));
PSFz.Ispline=single(getsmoothspline(Iz(mp-dd:mp+dd,mp-dd:mp+dd,:),ph));
PSFz.frequency=ph.frequency;
PSFz.phaseshifts=phaseshifts;
PSFz.factor=ones(4,1);
PSFz.factor([2 4])=factor;

validatemodel(PSFz,ph,'fitzast_<AB_i>') %z aligned by fitted z_astig


[imstackalignedpn,factor]=normalizequadrants(imstackalignedp);
[Im,Am,Bm,PSFm]=make4Pimodel(squeeze(mean(imstackalignedpn(:,:,:,:,:),4)),phaseshifts,ph.frequency,ph);

validatemodel(PSFm,ph,'fitzph_<PSF>') %z aligned by z_phase

%not better. Redo Transformation with fitted localizations?


%now testing needed. Which of the IAB produce all quadrants of calibration stack faithfully?
%use all for fitting, then look at residuals in quadrants
%look at phase etc.

%determine average x,y,z,phase for each bead for each channel. restrict to
%center?



    %cut out rois
    %fit unlinked
    %determine true x,y,z,phi for each bead


% find bead pairs
% cut out beads and move according to fit to perfect overlap
% calculate I, A, B for every bead
% perform registerPSF on I,A,B and all 4 channels together
% tilted coverslip: do we need to adjust phase for every bead to have same
% A, B, I? Use fitted phase phi for htis?

%alternatively: 
%as before for 2 channels average a 4 channel PSF (with fringes). Then use
%this for A,B,I



end


function [stack,filenumber,dT]=bead2stack(beads)
ss=size(beads(1).stack.image);
if length(ss)==3
    stack=zeros(ss(1),ss(2),ss(3),length(beads));
    for k=length(beads):-1:1
        stack(:,:,:,k)=beads(k).stack.image;
        filenumber(k)=beads.filenumber;
    end
elseif length(ss)==4
    stack=zeros(ss(1),ss(2),ss(3),length(beads),ss(4));
    numpar=6;
    dT=zeros(numpar,ss(4),ss(3),length(beads));
    for k=length(beads):-1:1
        stack(:,:,:,k,:)=beads(k).stack.image;
        filenumber(k)=beads.filenumber;
        for zz=1:ss(3)
            dT(1:2,:,zz,k)=squeeze(beads(k).shiftxy(1,[2 1],:));
        end
    end    
end
end

function [phaseshiftso,frequencyo]=getphaseshifts(allPSFs,ax,p)
ss=size(allPSFs);
range=(ss(1)+1)/2+[-1 1];
fw=20;
fw=ceil(500/p.dz);
frange=round(ss(3)/2+(-fw:fw)');

f=(1:ss(3))'-ss(3)/2;
intall=[];
%  kapprox=1*pi*4/max(f);
for k=1:ss(4)
    inth=squeeze(sum(sum(allPSFs(range,range,:,k),1),2));
    intall(:,k)=inth;
end
normn=sum(intall,2)/2;   % normalize by i1+i2+i3+i4
intnf=intall./normn;
intn=intnf(frange,:);
%k phi1 phi2 As Bs...
    lb1=[(max(intn(:))-min(intn(:)))/4 -inf -inf  min(intn(:)) -inf -inf ];
    ub1=[max(intn(:)) inf 0  max(intn(:)) inf 0 ];
    s1=[(max(intn(:))-min(intn(:)))/2 0 0 0.5 0 0];
    [~,indmax]=max(intn(:,1));
    
    
    %find kapprox
    indzero=find(f>=0,1,'first');
    inttest=intnf(indzero:end,1);
    ind1=find(inttest>=0.5,1,'first');
    inttest2=inttest(ind1:end);
    ind2=find(inttest2<=0.5,1,'first');
    inttest3=inttest2(ind2:end);
    ind3=find(inttest3>=0.5,1,'first');
     inttest4=inttest3(ind3:end);
    ind4=find(inttest4<=0.5,1,'first');   
    kapprox=pi/(ind3+ind4)*2
    
%     st1=[kapprox 0 0.5 0 0 0.5 0 0];
%     lba1=[0 -pi]
%     for k=1:size(intn,2)
%         fitpg=lsqcurvefit(@zintp,st1,f(frange),intn(:,k),lba1,uba1,[],0);
%     end
%     
%     
     phasestart1=pi/2-indmax*kapprox+pi; if phasestart1<0,phasestart1=phasestart1+2*pi;end;

lba=horzcat(-inf,-inf,-inf,lb1,lb1,lb1,lb1);
uba=horzcat(inf,inf,inf,ub1,ub1,ub1,ub1);

% phasestart1=0;
% kapprox=0.5;
startpa=[kapprox,phasestart1, phasestart1+pi/2,s1,s1,s1,s1];
fitAB=0;
fitpg=lsqcurvefit(@zintpg,startpa,f(frange),intn,lba,uba,[],fitAB);
% 
% fitted0=zintpg(fitpg,f(frange),fitAB);
% fitted0=reshape(fitted0,length(frange),4);

fitAB=1;
fitpg2=lsqcurvefit(@zintpg,fitpg,f(frange),intn,lba,uba,[],fitAB);



fitted=zintpg(fitpg2,f(frange),fitAB);
fitted=reshape(fitted,length(frange),4);

hold (ax,'off')
plot(ax,f,intnf,':+')
% plot(ax,f(frange),intn)
hold(ax,'on')
fst=zintpg(startpa,f(frange));
% plot(ax,f(frange),fst(:,:),'b--')

% plot(ax,f(frange)',fitted0','r');
plot(ax,f(frange)',fitted','k');

phaseshiftso=fitpg2([2 3]);
frequencyo=fitpg2(1)/2;
title(ax,['frequency: ' num2str(frequencyo,3) ', phaseshift/pi: ' num2str(mod((phaseshiftso(2)-phaseshiftso(1))/pi,2),3)])
%   fnc=@(k,phi1,phi2,A11,A21,A31,B11,B21,B31,A12,A22,A32,B12,B22,B32,A13,A23,A33,B13,B23,B33,A14,A24,A34,B14,B24,B34,x) zintg(k,phi1,A11,A21,A31,B11,B21,B31,phi2,A12,A22,A32,B12,B22,B32,phi3,A13,A23,A33,B13,B23,B33,phi4,A14,A24,A34,B14,B24,B34,x);

end

function into=zintpg(p,xdat,fitAB)
if nargin <3
    fitAB=1;
end
into=horzcat(zintp([p(1) p(2) p(4:9)],xdat,fitAB),zintp([p(1) p(3) p(10:15)],xdat,fitAB),zintp([p(1) p(2)+pi p(16:21)],xdat,fitAB),zintp([p(1) p(3)+pi p(22:27)],xdat,fitAB));
% into=vertcat(zint(x,k,phi1,A11,A21,A31,B11,B21,B31),zint(x,k,phi2,A12,A22,A32,B12,B22,B32),zint(x,k,phi3,A13,A23,A33,B13,B23,B33),zint(x,k,phi4,A14,A24,A34,B14,B24,B34));
end
function into=zint(f,k,phi,A1,A2,A3,B1,B2,B3,fitAB)
Bg=B1+B2*f*fitAB+B3*f.^2*fitAB;
Am=A1+A2*f*fitAB+A3*f.^2*fitAB;
os=sin(k*f+phi);
into=Bg+Am.*os;
end

function into=zintp(p,f,fitAB)
if nargin <3
    fitAB=1;
end
k=p(1);phi=p(2);A1=p(3);A2=p(4)*fitAB;A3=p(5)*fitAB;B1=p(6);B2=p(7)*fitAB;B3=p(8)*fitAB;
into=zint(f,k,phi,A1,A2,A3,B1,B2,B3,fitAB);
% Bg=B1+B2*f+B3*f.^2;
% Am=A1+A2*f+A3*f.^2;
% os=sin(k*f+phi);
% into=Bg+Am.*os;
end

% function [I,A,B,PSF]=make4Pimodel(allPSFs,phaseshifts,frequency,p)
% %re-weight every PSF by relative transmission?
%     I1=(allPSFs(:,:,:,1)+allPSFs(:,:,:,3))/2;
%     I2=(allPSFs(:,:,:,2)+allPSFs(:,:,:,4))/2;
%     Iall=(I1+I2)/2;
%     
%     z=(1:size(allPSFs,3))'-round(size(allPSFs,3)/2);
% [A12,B12]=makeAB(allPSFs(:,:,:,1),allPSFs(:,:,:,2),Iall,z,frequency,phaseshifts(1),phaseshifts(2));
% [A41,B41]=makeAB(allPSFs(:,:,:,4),allPSFs(:,:,:,1),Iall,z,frequency,phaseshifts(4),phaseshifts(1));
% [A23,B23]=makeAB(allPSFs(:,:,:,2),allPSFs(:,:,:,3),Iall,z,frequency,phaseshifts(2),phaseshifts(3));
% [A34,B34]=makeAB(allPSFs(:,:,:,3),allPSFs(:,:,:,4),Iall,z,frequency,phaseshifts(3),phaseshifts(4));
% A=(A12+A23+A34+A41)/4;
% B=(B12+B23+B34+B41)/4;
% I=Iall;
% 
% if nargin>3
% PSF=makeIABspline(I,A,B,p);
% PSF.frequency=frequency;
% PSF.phaseshifts=phaseshifts;
% end
% 
% end
function [out,normf]=makeIABspline(I,A,B,p)

% normalize to central frames of I in 5 x 5 region


mp=ceil((size(A,1)-1)/2);
dd=floor((p.ROIxy-1)/2);
intz=squeeze(sum(sum(I(mp-2:mp+2,mp-2:mp+2,:),1),2));
normf=max(intz);

out.Aspline=single(getsmoothspline(A(mp-dd:mp+dd,mp-dd:mp+dd,:)/normf,p));
out.Bspline=single(getsmoothspline(B(mp-dd:mp+dd,mp-dd:mp+dd,:)/normf,p));
out.Ispline=single(getsmoothspline(I(mp-dd:mp+dd,mp-dd:mp+dd,:)/normf,p));
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

function cspline=getsmoothspline(V,p)

if ~isfield(p,'pixelsize')
    pixelsizeh=100;
else
    pixelsizeh=p.pixelsize{1}(1);
end
lambdax=p.smoothxy/pixelsizeh/100000;
lambdaz=p.smoothz/p.dz*100;
lambda=[lambdax lambdax lambdaz];
b3_0t=bsarray(double(V),'lambda',lambda);
zhd=1:1:b3_0t.dataSize(3);
dxxhd=1;
[XX,YY,ZZ]=meshgrid(1:dxxhd:b3_0t.dataSize(1),1:dxxhd:b3_0t.dataSize(2),zhd);
corrPSFhdt = interp3_0(b3_0t,XX,YY,ZZ,0);
cspline = Spline3D_interp(corrPSFhdt);
end

function [imstackn,n]=normalizequadrants(imstack)
if length(size(imstack))==5
    imstackh=squeeze(mean(imstack,4));
else
    imstackh=imstack;
end
ry=5:size(imstackh,1)-4;
rx=5:size(imstackh,2)-4;
rz=7:size(imstackh,3)-6;

imstack1=imstackh(ry,rx,rz,1)+imstackh(ry,rx,rz,3);
imstack2=imstackh(ry,rx,rz,2)+imstackh(ry,rx,rz,4);
m1=max(imstack1(:));m2=max(imstack2(:));
indg=imstack1>m1/2&imstack2>m2/2;
% sum(imstack1(indg))/sum(imstack2(indg))
n=mean(imstack1(indg)./imstack2(indg),'omitnan');
imstackn=imstack;

imstn=(imstack1+imstack2)/2;
rzs=round(size(imstn,3)/2)+[-1 1];
amp=sum(sum(mean(mean(imstn(:,:,rzs,:),4),3),2),1);

if length(size(imstack))==5
    imstackn(:,:,:,:,[2 4])=imstackn(:,:,:,:,[2 4])*n;
else
    imstackn(:,:,:,[2 4])=imstackn(:,:,:,[2 4])*n;
end

imstackn=imstackn/amp;

end

function [allPSFs,shiftedstack,corrout]=PSFcorrelation(beads,ph)
%for each channel: calculate average PSF after alignment 
%this is approximate, as different beads can have different phases 
sstack=size(beads(1).stack.image);
xbeads=getFieldAsVectorInd(beads,'pos',1);
fw=round((ph.zcorrframes-1)/2);
framerange=round(sstack(3)/2-fw:sstack(3)/2+fw);
corrout.numchannels=length(ph.settings_3D.y4pi);
% allPSFs=zeros(sstack(1),sstack(2),sstack(3),corrout.numchannels);
% shiftedstack=[];
for k=1:corrout.numchannels
    %calculate average PSF
    w4pi=ph.settings_3D.width4pi;
    indbh=xbeads>=(k-1)*w4pi+1 & xbeads<= k*w4pi;
    [allstacks,corrout.filenumber]=bead2stack(beads(indbh));
    [~,shiftedstackh,shift,indgood]=registerPSF3D_g(allstacks,[],struct('framerange',framerange,'normalize',false));
    % determine true position of the beads in the four channels
    xposh=getFieldAsVectorInd(beads(indbh),'pos',1)';
    yposh=getFieldAsVectorInd(beads(indbh),'pos',2)';
    
    beadtrue{k}(:,1)=xposh(indgood)-shift(indgood,2)'; %experimentally: this works :)
    beadtrue{k}(:,2)=yposh(indgood)-shift(indgood,1)';
%     corrout.beadtrue{k}(:,1)=xposh(indgood)-shift(indgood,2)'; %experimentally: this works :)
%     corrout.beadtrue{k}(:,2)=yposh(indgood)-shift(indgood,1)';
    %beads in all stacks could be shifted to different heights. Try to
    %compensate by subtracting mean shift~
    beadtrue{k}(:,3)=shift(indgood,3)';%-mean(shift(indgood,3)); %this is not yet tested, could be minus  
%     corrout.zshift0(k)=mean(shift(indgood,3));
    shiftxy{k}=shift(indgood,1:2);
    shiftedstack{k}=shiftedstackh(:,:,:,indgood);
end
ph.plottransform=false;
transformation=make4PiTransform(beadtrue,ph);
iAtot=(1:size(beadtrue{1},1))';
for k=2:corrout.numchannels
    post=transformation.transformToReference(k,beadtrue{k});
    [iA{k},iB{k}]=matchlocs(beadtrue{1}(:,1),beadtrue{1}(:,2),post(:,1),post(:,2),[0 0],3);
    iAtot=intersect(iAtot,iA{k});
end

indf{1}=iAtot;
for k=2:corrout.numchannels
    [~,ia,ib]=intersect(iAtot, iA{k});
    indf{k}=iB{k}(ib);
end

allPSFs=zeros(sstack(1),sstack(2),sstack(3),corrout.numchannels);
for k=1:corrout.numchannels
    corrout.beadtrue{k}=beadtrue{k}(indf{k},:);
    shiftedstack{k}=shiftedstack{k}(:,:,:,indf{k});
    allPSFs(:,:,:,k)=nanmean(shiftedstack{k},4);
    corrout.zshift0(k)=nanmean(corrout.beadtrue{k}(:,3));
    corrout.xyshift0(k,1:2)=squeeze(nanmean(shiftxy{k}(indf{k},:),1)) ;
end
% corrout.shiftedstack=shiftedstack;
end


% function [allPSFs,shiftedstack,corrout]=PSFcorrelation(beads,ph)
% %for each channel: calculate average PSF after alignment 
% %this is approximate, as different beads can have different phases 
% sstack=size(beads(1).stack.image);
% xbeads=getFieldAsVectorInd(beads,'pos',1);
% fw=round((ph.zcorrframes-1)/2);
% framerange=round(sstack(3)/2-fw:sstack(3)/2+fw);
% corrout.numchannels=length(ph.settings_3D.y4pi);
% allPSFs=zeros(sstack(1),sstack(2),sstack(3),corrout.numchannels);
% % shiftedstack=[];
% for k=1:corrout.numchannels
%     %calculate average PSF
%     w4pi=ph.settings_3D.width4pi;
%     indbh=xbeads>=(k-1)*w4pi+1 & xbeads<= k*w4pi;
%     [allstacks,corrout.filenumber]=bead2stack(beads(indbh));
%     [allPSFs(:,:,:,k),shiftedstackh,shift,indgood]=registerPSF3D_g(allstacks,[],struct('framerange',framerange,'normalize',ph.normalize));
%     % determine true position of the beads in the four channels
%     xposh=getFieldAsVectorInd(beads(indbh),'pos',1)';
%     yposh=getFieldAsVectorInd(beads(indbh),'pos',2)';
%     corrout.beadtrue{k}(:,1)=xposh(indgood)-shift(indgood,2)'; %experimentally: this works :)
%     corrout.beadtrue{k}(:,2)=yposh(indgood)-shift(indgood,1)';
%     %beads in all stacks could be shifted to different heights. Try to
%     %compensate by subtracting mean shift
%     corrout.beadtrue{k}(:,3)=shift(indgood,3)';%-mean(shift(indgood,3)); %this is not yet tested, could be minus  
%     corrout.zshift0(k)=mean(shift(indgood,3));
%     corrout.xyshift0(k,1:2)=squeeze(mean(shift(indgood,1:2),1));
%     shiftedstack{k}=shiftedstackh(:,:,:,indgood);
% end
% 
% 
% 
% end

function transform=make4PiTransform(beadtrue,ph)
%calculate transformN
transform=interfaces.LocTransformN;
pt.mirror=0; %mirror already taken care of when reading in images
settings_3D=ph.settings_3D;
pt.xrange=[1 settings_3D.width4pi];
pt.yrange=[1 settings_3D.height4pi];
pt.unit='pixel';
pt.type='projective';
transform.setTransform(1,pt)
iAaa=1:size(beadtrue{1},1);
if ~isfield(ph,'plottransform') || ph.plottransform
th=uitab(ph.tabgroup,'Title','transform');
tabgroup=uitabgroup(th);
ploton=true;
else 
    ploton=false;
end
for k=2:length(beadtrue)
    pt.xrange=[(k-1)*settings_3D.width4pi+1 k*settings_3D.width4pi];
    transform.setTransform(k,pt)
    if ploton
        tab=(uitab(tabgroup,'Title',['T' num2str(k)]));ph.ax=axes(tab);
    else
        ph.ax=[];
    end
    [~ ,iAa,iBa]=transform_locs_simpleN(transform,1, beadtrue{1},k,beadtrue{k},ph); 
    %extend transform locs by iterative transform - remove outliers. As
    %done for normal calibrator.
    iAaa=intersect(iAa,iAaa);
end
end



function [img,beads]=validatemodel(PSF,ph,titlet)
if nargin<3
    titlet='results';
end
Nfree=false;
ph.isglobalfit=true;

[beads,ph]=images2beads_globalfitN(ph); 
   %%%%%%%%%%%%%%%%%
   beads=beads(1)
%%%%%%%%%%
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
imstacksq=imsqueeze(ph.rangeh, ph.rangeh, :, :);
iterations=50;
z0=ph.zstart;
[P,CRLB1 LL] = mleFit_LM_4Pi(single(imstacksq(:, :, :, :)),uint32(shared),iterations,single(PSF.Ispline), single(PSF.Aspline),single(PSF.Bspline),single(dTAll),single(ph.phi0),z0);

% [P,CRLB1 LL] = CPUmleFit_LM_MultiChannel_4pi(single(imstacksq(:, :, :, :)),uint32(shared),iterations,single(PSF.Ispline), single(PSF.Aspline),single(PSF.Bspline),single(dTAll),single(ph.phi0),z0);

img.imstacksq=imstacksq;
img.sim=sim;
img.fit.P=P;
img.fit.CRLB=CRLB1;
img.fit.PSF=PSF;
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

% calculate residuals for each bead. Should it beased on average x,y,z in center?

end