% challengescript
% parameters
%% DHNPC 10nm  C1_MT1.N1.LD\DHNPC

% 44.453 45.166 9.744
if 0
    dx=44.453; %corrections from bead fit.
    dy=45.166;
    dz=9.744;
    photonfactor=0.662;
    
    %% filter used
    
    groupfilter = 1;
    LLfilter = 1;
    iterationfilter = 1;
    filterint=1;
    boarderfilter = 1;
    clusterfilter = 1;
    zfilter = 0;
    
    
    %density
    densitysize_xy=50;
    densitysize_z=100;
    densitycutoff=20;
    zmin=-800;zmax=800;
    bgmin=0; bgmax=14000; %background filter
    locprec_cutoff=60;
    locprecz_cutoff=100;
    phot_cutoff=0;
    LLrel_cutoff=-2.2;
    group_dT=0;
    border=5; %distance from min/max: if fit did converge to border
end



%% DHNPC10nm c3_MT3.N2.LD\DHNPC

% 44.453 45.166 9.744
if 0
    dx=44.453; %corrections from bead fit.
    dy=45.166;
    dz=9.744;
    photonfactor=0.662;
    
    %% filter used
    
    groupfilter = 1;
    LLfilter = 1;
    iterationfilter = 1;
    filterint=1;
    boarderfilter = 1;
    clusterfilter = 1;
    zfilter = 0;
    
    %density
    densitysize_xy=50;
    densitysize_z=100;
    densitycutoff=10;
    zmin=-800;zmax=800;
    bgmin=0; bgmax=11; %background filter
    locprec_cutoff=80;
    locprecz_cutoff=140;
    phot_cutoff=.500;
    LLrel_cutoff=-1.5;
    group_dT=0;
    border=5; %distance from min/max: if fit did converge to border
end


%% DHNPC 10nm  C2_MT2.N1.HD\DHNPC

% 44.453 45.166 9.744
if 0
    dx=44.453; %corrections from bead fit.
    dy=45.166;
    dz=9.744;
    photonfactor=0.662;
    
    %% filter used
    
    groupfilter = 1;
    LLfilter = 1;
    iterationfilter = 1;
    filterint=1;
    boarderfilter = 1;
    clusterfilter = 1;
    zfilter = 0;
    
    
    %density
    densitysize_xy=50;
    densitysize_z=250;
    densitycutoff=15;
    zmin=-800;zmax=800;
    bgmin=0; bgmax=14000; %background filter
    locprec_cutoff=30;
    locprecz_cutoff=60;
    phot_cutoff=0;
    LLrel_cutoff=-7;
    group_dT=0;
    border=5; %distance from min/max: if fit did converge to border
end

%% DHNPC 10nm  \C4_MT4.N2.HD\DHNPC

% 44.453 45.166 9.744
if 0
    dx=44.453; %corrections from bead fit.
    dy=45.166;
    dz=9.744;
    photonfactor=0.662;
    
    %% filter used
    
    groupfilter = 1;
    LLfilter = 1;
    iterationfilter = 1;
    filterint=1;
    boarderfilter = 1;
    clusterfilter = 1;
    zfilter = 0;
    
    
    %density
    densitysize_xy=50;
    densitysize_z=250;
    densitycutoff=9;
    zmin=-800;zmax=800;
    bgmin=0; bgmax=14000; %background filter
    locprec_cutoff=80;
    locprecz_cutoff=150;
    phot_cutoff=0;
    LLrel_cutoff=-1.8;
    group_dT=0;
    border=5; %distance from min/max: if fit did converge to border
end



%% BP 10nm C1_MT1.N1.LD

% -54.864 45.095 12.802

if 0
    dx=-54.864; %corrections from bead fit.
    dy=45.095;
    dz=12.802;
    photonfactor=1.6732;
    
    %% filter used
    
    groupfilter = 1;
    LLfilter = 1;
    iterationfilter = 1;
    filterint=1;
    boarderfilter = 1;
    clusterfilter =1;
    zfilter = 0;
    
    %density
    densitysize_xy=50;
    densitysize_z=50;
    densitycutoff=3;
    zmin=-800;zmax=800;
    bgmin=75; bgmax=125; %background filter
    locprec_cutoff=40;
    locprecz_cutoff=80;
    phot_cutoff=100;
    LLrel_cutoff=-3;
    group_dT=0;
    border=5; %distance from min/max: if fit did converge to border
end


%% BP 10nm c3_MT3.N2.LD

% -54.864 45.095 12.802

if 0
    dx=-54.864; %corrections from bead fit.
    dy=45.095;
    dz=12.802;
    photonfactor=1.6732;
    
    %% filter used
    
%     groupfilter = 0;
%     LLfilter = 0;
%     iterationfilter = 1;
%     filterint=1;
%     boarderfilter = 1;
%     clusterfilter =1;
%     zfilter = 0;
%   
%     
%     %density
%     densitysize_xy=50;
%     densitysize_z=200;
%     densitycutoff=10;
%     zmin=-800;zmax=800;
%     bgmin=7; bgmax=10; %background filter
%     locprec_cutoff=40;
%     locprecz_cutoff=100;
%     phot_cutoff=150;
%     LLrel_cutoff=-4.5;
%     group_dT=0;
%     border=5; %distance from min/max: if fit did converge to border

    groupfilter = 1;
    LLfilter = 1;
    iterationfilter = 1;
    filterint=1;
    boarderfilter = 1;
    clusterfilter =1;
    zfilter = 0;
    
    %density
    densitysize_xy=50;
    densitysize_z=50;
    densitycutoff=2;
    zmin=-800;zmax=800;
    bgmin=6.5; bgmax=10.5; %background filter
    locprec_cutoff=60;
    locprecz_cutoff=120;
    phot_cutoff=.100;
    LLrel_cutoff=-2.6;
    group_dT=0;
    border=5; %distance from min/max: if fit did converge to border
end

%% BP 10nm C2_MT2.N1.HD\BP

% -54.864 45.095 12.802

if 0
    dx=-54.864; %corrections from bead fit.
    dy=45.095;
    dz=12.802;
    photonfactor=1.6732;
    
    %% filter used
    
    groupfilter = 1;
    LLfilter = 1;
    iterationfilter = 1;
    filterint=1;
    boarderfilter = 1;
    clusterfilter =1;
    zfilter = 0;
    
    %density
    densitysize_xy=50;
    densitysize_z=50;
    densitycutoff=3;
    zmin=-800;zmax=800;
    bgmin=.75; bgmax=12115; %background filter
    locprec_cutoff=60;
    locprecz_cutoff=100;
    phot_cutoff=.100;
    LLrel_cutoff=-10;
    group_dT=0;
    border=5; %distance from min/max: if fit did converge to border
end

%% BP 10nm \C4_MT4.N2.HD\BP

% -54.864 45.095 12.802

if 0
    dx=-54.864; %corrections from bead fit.
    dy=45.095;
    dz=12.802;
    photonfactor=1.6732;
    
    %% filter used
    
    groupfilter = 1;
    LLfilter = 1;
    iterationfilter = 1;
    filterint=1;
    boarderfilter = 1;
    clusterfilter =1;
    zfilter = 0;
    
    %density
    densitysize_xy=50;
    densitysize_z=50;
    densitycutoff=2;
    zmin=-800;zmax=800;
    bgmin=.75; bgmax=12115; %background filter
    locprec_cutoff=80;
    locprecz_cutoff=160;
    phot_cutoff=.100;
    LLrel_cutoff=-3.5;
    group_dT=0;
    border=5; %distance from min/max: if fit did converge to border
end


%% AS 10nm C1_MT1.N1.LD
if 0
    dx=-4.604; %corrections from bead fit.
    dy=95.165;
    dz=3.150;
    photonfactor=0.872;
    
    %% filter used
    
    groupfilter = 1;
    LLfilter = 1;
    iterationfilter = 1;
    filterint=1;
    boarderfilter = 1;
    clusterfilter = 1;
    zfilter = 0;
    
    
    %density
    densitysize_xy=50;% for clusterfilter
    densitysize_z=100;% for clusterfilter
    densitycutoff=3;% for clusterfilter
    zmin=-800;zmax=800;% for zfilter
    bgmin=-inf; bgmax=inf; %background filter
    locprec_cutoff=60; %groupfilter
    locprecz_cutoff=120; %groupfilter
    phot_cutoff=.100; %groupfilter
    LLrel_cutoff=-2; %LLfilter
    group_dT=0; %regroup
    border=5; %distance from min/max: if fit did converge to border; boarderfilter
end

%% AS 10nm C3_MT3.N1.LD
if 0
    dx=-4.604; %corrections from bead fit.
    dy=95.165;
    dz=3.150;
    photonfactor=0.872;
    
    %% filter used
    
    groupfilter = 1;
    LLfilter = 1;
    iterationfilter = 1;
    filterint=1;
    boarderfilter = 1;
    clusterfilter = 1;
    zfilter = 0;
    
    
    %density
    densitysize_xy=80;% for clusterfilter
    densitysize_z=100;% for clusterfilter
    densitycutoff=3;% for clusterfilter
    zmin=-800;zmax=800;% for zfilter
    bgmin=-inf; bgmax=inf; %background filter
    locprec_cutoff=80; %groupfilter
    locprecz_cutoff=160; %groupfilter
    phot_cutoff=.100; %groupfilter
    LLrel_cutoff=-2; %LLfilter
    group_dT=0; %regroup
    border=5; %distance from min/max: if fit did converge to border; boarderfilter
end


%% AS 10nm C2_MT2.N1.HD_version2_sig1.5_roi5_
if 0
    dx=-4.604; %corrections from bead fit.
    dy=95.165;
    dz=3.150;
    photonfactor=0.872;
    
    %% filter used
    
    groupfilter = 1;
    LLfilter = 1;
    iterationfilter = 1;
    filterint=0;
    boarderfilter = 0;
    clusterfilter = 1;
    zfilter = 0;
    
    
    %density
    densitysize_xy=70;
    densitysize_z=130;
    densitycutoff=20;
    zmin=-800;zmax=800;
    bgmin=65; bgmax=175; %background filter
    locprec_cutoff=50;
    locprecz_cutoff=100;
    phot_cutoff=100;
    LLrel_cutoff=-1.8;
    group_dT=0;
    border=5; %distance from min/max: if fit did converge to border
end

%% AS 10nm C4_MT4.N2.HD_version2_sig1_roi5_
if 0
    dx=-4.604; %corrections from bead fit.
    dy=95.165;
    dz=3.150;
    photonfactor=0.872;
    
    %% filter used
    
    groupfilter = 1;
    LLfilter = 1;
    iterationfilter = 1;
    filterint=1;
    boarderfilter = 1;
    clusterfilter = 1;
    zfilter = 0;
    
    
    %density
    densitysize_xy=80;% for clusterfilter
    densitysize_z=160;% for clusterfilter
    densitycutoff=2;% for clusterfilter
    zmin=-800;zmax=800;% for zfilter
    bgmin=-inf; bgmax=inf; %background filter
    locprec_cutoff=80; %groupfilter
    locprecz_cutoff=160; %groupfilter
    phot_cutoff=.100; %groupfilter
    LLrel_cutoff=-2; %LLfilter
    group_dT=0; %regroup
    border=5; %distance from min/max: if fit did converge to border; boarderfilter
end




%% 2D C5_ER1.N3.LD Smooth4
if 0
    dx=7.688; %corrections from bead fit.
    dy=-15.484;
    dz=-0.443+25;
    photonfactor=0.896;
    
    %% filter used
    
    groupfilter = 1;
    LLfilter = 1;
    iterationfilter = 1;
    filterint=1;
    boarderfilter = 1;
    clusterfilter = 0;
    zfilter = 0;
    
    
    %density
    densitysize_xy=50;
    densitysize_z=1000;
    densitycutoff=2;
    zmin=-800;zmax=800;
    bgmin=-inf; bgmax=inf; %background filter
    locprec_cutoff=80;
    locprecz_cutoff=160;
    phot_cutoff=0.1;
    LLrel_cutoff=-2.5;
    group_dT=0;
    border=5; %distance from min/max: if fit did converge to border
end

%% 2D \c3_MT3.N2.LD\2D Smooth4
if 0
    dx=7.688; %corrections from bead fit.
    dy=-15.484;
    dz=-0.443+25;
    photonfactor=0.896;
    
    %% filter used
    
    groupfilter = 1;
    LLfilter = 1;
    iterationfilter = 1;
    filterint=1;
    boarderfilter = 1;
    clusterfilter = 0;
    zfilter = 0;
    
    
    %density
    densitysize_xy=100;
    densitysize_z=300;
    densitycutoff=2;
    zmin=-800;zmax=800;
    bgmin=-inf; bgmax=inf; %background filter
    locprec_cutoff=80;
    locprecz_cutoff=160;
    phot_cutoff=0.1;
    LLrel_cutoff=-1.4;
    group_dT=0;
    border=5; %distance from min/max: if fit did converge to border
end

%% 2D \C4_MT4.N2.HD\2D Smooth4
if 1
    dx=7.688; %corrections from bead fit.
    dy=-15.484;
    dz=-0.443+25;
    photonfactor=0.896;
    
    %% filter used
    
    groupfilter = 0;
    LLfilter = 1;
    iterationfilter = 1;
    filterint=1;
    boarderfilter = 0;
    clusterfilter = 1;
    zfilter = 0;
    
    
    %density
    densitysize_xy=100;
    densitysize_z=1000;
    densitycutoff=100;
    zmin=-800;zmax=800;
    bgmin=-inf; bgmax=inf; %background filter
    locprec_cutoff=100;
    locprecz_cutoff=400;
    phot_cutoff=0.1;
    LLrel_cutoff=-2;
    group_dT=0;
    border=5; %distance from min/max: if fit did converge to border
end

%% 2D \C6_ER2.N3.HD\2D Smooth4
if 1
    dx=7.688; %corrections from bead fit.
    dy=-15.484;
    dz=-0.443+25;
    photonfactor=0.896;
    
    %% filter used
    
    groupfilter = 1;
    LLfilter = 1;
    iterationfilter = 1;
    filterint=1;
    boarderfilter = 0;
    clusterfilter = 1;
    zfilter = 0;
    
    
    %density
    densitysize_xy=100;
    densitysize_z=1000;
    densitycutoff=50;
    zmin=-800;zmax=800;
    bgmin=-inf; bgmax=inf; %background filter
    locprec_cutoff=100;
    locprecz_cutoff=500;
    phot_cutoff=0.1;
    LLrel_cutoff=-2;
    group_dT=0;
    border=5; %distance from min/max: if fit did converge to border
end



%%
compare=true; %open compare java
% %% fit time
% fn=g.locData.files.file.name;
% l=load(fn);
% fitt=l.saveloc.fitparameters.processfittime;
% disp(['fit time without loading: ' num2str(fitt,3) ' s']);
%%

ld=g.locData;
maxiter=max(ld.loc.iterations);
% group_dX=[15 30 50] 


PSF_cutoff_min=75; %only for Gaussian, so not relevant?
PSF_cutoff_max=200;


filenumber =1;

[~,filename]=fileparts(ld.files.file(filenumber).name);


group_dX=min(50,2*median(ld.loc.locprecnm));

% determine HD, LD, bright, dark
%regroup with smaller dx, dt (depending on brightness, density)
ld.regroup(group_dX,group_dT)

%filter
indgood=true(size(ld.loc.xnm));
%  combined PSF grouping filter. PSF not relevant?

if groupfilter
    indgroupfilter=((...
        ld.loc.locprecnm<locprec_cutoff...
        & ld.loc.locprecznm<locprecz_cutoff...
        & ld.loc.PSFxnm>PSF_cutoff_min & ld.loc.PSFxnm<PSF_cutoff_max ...
        & ld.loc.bg>bgmin & ld.loc.bg<bgmax  ...
        & ld.loc.phot>phot_cutoff ...
        )| ld.loc.numberInGroup>1);
    indgood=indgood & indgroupfilter;
    disp(['groupfilter: ' num2str(sum(~indgroupfilter)/length(indgroupfilter)*100) ' rem: ' num2str(sum(indgood))])
end


%  LL filtering
if LLfilter
    if isfield(ld.loc,'LLrel') %do this filtr before regrouping?
        indll=ld.loc.LLrel>LLrel_cutoff;
    end
    outnold=sum(~indgood);
    indgood=indgood & indll;
    out=sum(~indgood)-outnold;
    disp(['LLrel: ' num2str(out/length(indgood)*100) ' rem: ' num2str(sum(indgood))])
end
% disp(['groupfilter: ' num2str(1-sum(indgroupfilter)/length(indgroupfilter))])
%iterations
if iterationfilter
    inditer=ld.loc.iterations<maxiter;
    
    outnold=sum(~indgood);
    indgood=indgood & inditer;
    out=sum(~indgood)-outnold;
    disp(['iter: ' num2str(out/length(indgood)*100) ' rem: ' num2str(sum(indgood))])
end


if filterint
%stripe artifacts: can be reduced if integer numbers are removed
indint=round(ld.loc.xnm)-ld.loc.xnm == 0 | ...
    round(ld.loc.ynm)-ld.loc.ynm == 0 | ...
    round(ld.loc.znm)-ld.loc.znm == 0;
indround= ~indint;


outnold=sum(~indgood);
indgood=indgood & indround;
out=sum(~indgood)-outnold;
disp(['round: ' num2str(out/length(indgood)*100) ' rem: ' num2str(sum(indgood))])
end
%border filtering

if boarderfilter
indborder=ld.loc.xnm>min(ld.loc.xnm)+border & ...
    ld.loc.xnm<max(ld.loc.xnm)-border & ...
    ld.loc.ynm>min(ld.loc.ynm)+border & ...
    ld.loc.ynm<max(ld.loc.ynm)-border & ...
    ld.loc.znm>min(ld.loc.znm)+border & ...
    ld.loc.znm<max(ld.loc.znm)-border;
% indborder=indborder;

outnold=sum(~indgood);
% indgood=indgood & indborder;
out=sum(~indgood)-outnold;
disp(['border: ' num2str(out/length(indgood)*100) ' rem: ' num2str(sum(indgood))])
end
% only current filenumber:
indgood =indgood &ld.loc.filenumber==filenumber;

%idea: modified cluster density: dx,y,z depends on median local locprecnm
%  density calculation, removal of single localisations
fdcal=figure(233);
dcal=plugin('Analyze','cluster','density_calculator',fdcal,g.P);
dcal.attachLocData(ld);
dcal.makeGui;
p=dcal.getGuiParameters;
p.countwhat.Value=1;
p.countingsize_xy=densitysize_xy;
p.countingsize_z=densitysize_z;
dcal.setGuiParameters(p);
dcal.useind=indgood;
dcal.processgo;
indcluster=ld.loc.clusterdensity>=densitycutoff;

if clusterfilter
    outnold=sum(~indgood);
    indgood=indgood & indcluster;
    out=sum(~indgood)-outnold;
    disp(['cluster: ' num2str(out/length(indgood)*100) ' rem: ' num2str(sum(indgood))])
end


if zfilter
    %  edges: z min, z max: remove or leave?
    if isfield(ld.loc,'znm')
        indz = (ld.loc.znm>zmin & ld.loc.znm < zmax);
    end
    outnold=sum(~indgood);
    indgood=indgood & indz;
    out=sum(~indgood)-outnold;
    disp(['z: ' num2str(out/length(indgroupfilter))])
end

g.locData.setloc('challengefiltered',single(indgood));
g.locData.regroup(group_dX,group_dT);

% write x,y,z from grouped to ungrouped
copygroupfields={'xnm','ynm','znm'};
ldc=ld.copy; %dont overwrite in SMaP
ldc.removelocs(~indgood)
ldc.regroup(group_dX,group_dT);

if 1
[gi,sorti]=sort(ldc.loc.groupindex);
[gig,sortg]=sort(ldc.grouploc.groupindex);
inds=1;indg=1;
for k=1:gi(end)
    while gi(inds)<k
        inds=inds+1;
    end
    while gig(indg)<k
        indg=indg+1;
    end
    inds2=inds;
    while gi(inds2)==k && inds2<length(gi)
        inds2=inds2+1;
    end
    ind2g=indg;
    while gig(ind2g)==k && inds2<length(gig)
        ind2g=ind2g+1;
    end    
    indc=inds:inds2-1;
    gi(indc);
    gig(indg);
    if ~isempty(indc)
        if any(gi(indc)~=gig(indg))
            display('inconsitency group index')
        end
        for l=1:length(copygroupfields)
            if isfield(ldc.loc,copygroupfields{l})
                ldc.loc.(copygroupfields{l})(sorti(indc))=ldc.grouploc.(copygroupfields{l})(sortg(indg));
            end
        end
    end
end
end

if compare
fchallege=figure(234);
cs=plugin('Analyze','other','CompareToGroundTruthChallenge',fchallege,g.P);
cs.attachLocData(ldc);
cs.makeGui;
p=cs.getGuiParameters;
p.onlyfiltered=0;
p.offsetxyz=[dx dy dz];
p.photonfactor=photonfactor;
p.shiftframe=0;
cs.setGuiParameters(p);
cs.processgo;
end
disp('done')

% write  dx, dy, dz from bead fit into challenge compare plugin

%% create best J vs RMS plot
% GTfile='/Volumes/t2ries/projects/SMLMChallenge2018/T1_MT0.N1.LD/activations.csv';
% % GTfile='/Volumes/t2ries/projects/SMLMChallenge2018/T2_MT0.N1.LD/activations.csv';
% GTfile='Z:\projects\SMLMChallenge2018\T1_MT0.N1.LD\activations.csv';
% dat = csvread(GTfile,1,0);
% excessnoise=2;
% phot=dat(:,6)*.9/excessnoise;
% bg=90; 
% 
% PSF0=100; a=100; 
% 
% %PSF(z)
% lambda=600;
% w0=2*PSF0;
% z=dat(:,5);
% zR=pi*w0^2/lambda;
% wz=w0*sqrt(1+(z/zR).^2);
% PSF=wz/2;
% 
% locprecnm=sqrt((PSF.^2+a^2/12)./phot.*(16/9+8*pi*(PSF.^2+a^2/12)*bg./phot/a^2));
% lps=sort(locprecnm);
% norm=(1:length(lps))';
% lpsj=sqrt(cumsum(lps.^2)./norm);
% figure(88);
% subplot(1,2,2)
% hold off
% plot(norm/max(norm),lpsj,'.')
% hold on
% plot([0 1],[1 1]*lpsj(1))
% ax=gca;
% ax.YLim(1)=0;
% ax.YLim(2)=quantile(lpsj,0.99);
% xlabel('Recall')
% ylabel('RMS (nm)')
% title('best possible RMS vs recall')
% 
% subplot(1,2,1);
% histogram(locprecnm)
% xlabel('localization precision nm')
% xlim([0 quantile(locprecnm,0.98)])
% title('localization precision Mortenson')
