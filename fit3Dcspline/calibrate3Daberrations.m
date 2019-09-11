%  Copyright (c)2017 Ries Lab, European Molecular Biology Laboratory,
%  Heidelberg.
%  
% 
%  Additional permission under GNU GPL version 3 section 7
%  
%  If you modify this Program, or any covered work, by
%  linking or combining it with libraries required for interaction
%  with analysis programs such as Igor Pro or Matlab,
%  the licensors of this Program grant you additional permission
%  to convey the resulting work.
%%
function zcorr= calibrate3Daberrations(locs,pin)

%loc.filenumber
%loc.frame
%loc.phot
%loc.x
%loc.y
%loc.z

% default parameters. Those get overwritten by parameters in pin
p.glassframe=[]; %position of the frame in which beads on the glass are in focus. empty [] for automatic detection
p.dz=10;  %distance of objective positions in nm (in objective space, i.e. no refractive index mismatch correction)
p.smoothz=1/1000; %parameter how much to smooth the interpolation in z (fitted z direction).
p.smoothframe=2/10/p.dz; %smoothing factor along the frame direction
p.cutoffrefine=100; %maximum allowed distance of a beads z position from interpolated correction curve.
p.maxrange=800; %range around zero in which the correction is calcualted.

p=copyfields(p,pin); %overwrite with passed on parameters

f=figure('Name','Calibrate depth-induced aberrations');
p.tabgroup=uitabgroup(f);

%get beads from localizations
beads=segmentb_so(locs,p.dz); %find bead positions

% get true positions in unists of frames f0 for beads
for k=length(beads):-1:1
    [beads(k).f0]=getf0Z_so(beads(k).loc,p.dz);
end

%remove beads that resulted in none or wrong f0
f0all=([beads(:).f0]);
badind=f0all<min(locs.frame)|isnan(f0all); 
beads(badind)=[];

% determine position of the glass: no beads below glass
if isempty(p.glassframe)
    p.axhere=axes(uitab(p.tabgroup,'Title','f0'));
    p.glassframe=getf0glass(beads,p);
end
%correct all frame values by glassframe
for k=1:length(beads)
    beads(k).f0=beads(k).f0-p.glassframe;
    beads(k).loc.frame=beads(k).loc.frame-p.glassframe;
end

%calculate relevant other coordinates
axh=axes(uitab(p.tabgroup,'Title','focal plane vs fitted z'));
hold off
for k=1:length(beads)
    beads(k).loc.zobjective=(beads(k).loc.frame)*p.dz; %global normalized z-coordinate: objective position in nm above glass. old: .zglass. rename to zobjective
    beads(k).loc.zobjectiverelative=(beads(k).f0-beads(k).loc.frame)*p.dz; %z coordinate relative to true z-postion of bead. old: z0relative.%zobjective moved opposite to the relative bead position. Thus the minus in front of .frame
    beads(k).loc.dzcorr=beads(k).loc.zobjectiverelative-beads(k).loc.z;
    indplot=abs(beads(k).loc.zobjectiverelative)<p.maxrange;
    plot(axh,beads(k).loc.z(indplot),beads(k).loc.zobjective(indplot),'-')
    hold on
end
ylabel('focal plane')
xlabel('fitted z position (nm)')

p.axhere=[];
phere=p;
phere.smoothing=[0.05 0.002];
% iteratively determine spline approximation and remove beads that
% are too far away
[ZcorrInterp]=getZinterp(beads,[],phere);
phere=p;
phere.cutoffrefine=500;
[ZcorrInterp]=getZinterp(beads,ZcorrInterp,phere);
err1=geterrors(beads,ZcorrInterp);
goodind=find(true(length(beads),1));
beads2=beads;
while  1% length(beads2)>length(beads)/2
    cutoff=3*mean(err1,'omitnan');
    badind=(err1>cutoff|isnan(err1));
    if sum(badind)==0
        break
    end
    goodind=goodind(~badind);
    beads2=beads(goodind);
    [ZcorrInterp]=getZinterp(beads2,ZcorrInterp,p);
    %calculate errors
    err1=geterrors(beads2,ZcorrInterp);    
end

%take out beads for testing

numtestbeads=15;
f0=[beads2(:).f0];
f0test=linspace(min(f0),max(f0),numtestbeads)
for k=1:length(f0test)
    [~,testind(k)]=min(abs(f0test(k)-f0));

end
% numtotal=length(beads2);
% testind=round(numtotal*rand(numtestbeads,1));
beadtest=beads2(testind);
beads2(testind)=[];
%plot output
p.axhere=axes(uitab(p.tabgroup,'Title','Interpolation'));
[ZcorrInterp]=getZinterp(beads2,ZcorrInterp,p);
zcorr=ZcorrInterp.interp;  

ax2=axes(uitab(p.tabgroup,'Title','dz vs Zfit'));
zfit=-p.maxrange:p.maxrange;
minzobj=-p.maxrange;
objectivepos=linspace(minzobj,p.axhere.XLim(2),25);

col=jet(length(objectivepos));
for k=1:length(objectivepos)
    dzh=zcorr(objectivepos(k)*ones(size(zfit)),zfit);
    plot(ax2,zfit, (dzh),'Color',col(k,:))
    hold(ax2,'on');
%     plot(ax3,zfit, (dzh+zfit)./zfit,'Color',col(k,:))
%     hold(ax3,'on');
end
xlabel(ax2,'fitted z position (nm)');
ylabel(ax2,'correction dz (nm)');

%Validate: correct beads for testing
ax1=axes(uitab(p.tabgroup,'Title','Validation'));
f=ax1.Parent;
ax2=axes(f,'Position',[0.5 0 1 1]);
subplot(1,2,1,ax1);
subplot(1,2,2,ax2);

beads=beadtest;
for k=length(beads):-1:1
   dZ=zcorr(beads(k).loc.zobjective,beads(k).loc.z);
   beads(k).loc.zcorrected=beads(k).loc.z+dZ;
%    if ~any(goodind==k)
%        col='r.';
%        plot(ax1,beads(k).loc.zobjectiverelative,beads(k).loc.z,col)
%        plot(ax2,beads(k).loc.zobjectiverelative,beads(k).loc.zcorrected,col)
%        hold(ax1,'on');
%        hold(ax2,'on');
%    end
end
length(beads2)
color=jet(length(beads));
for k=length(beads):-1:1
    if any(goodind==k)
       goodz=abs(beads(k).loc.zobjectiverelative)<1000;
       plot(ax1,-1.*beads(k).loc.zobjectiverelative(goodz),-1.*beads(k).loc.z(goodz),'.','Color',color(k,:))
       plot(ax2,-1.*beads(k).loc.zobjectiverelative(goodz),-1.*beads(k).loc.zcorrected(goodz),'.','Color',color(k,:))
      
       hold(ax1,'on');
       hold(ax2,'on');
 
    end
end
plot(ax1,[-1000 1000],[-1000 1000],'k')
plot(ax2,[-1000 1000],[-1000 1000],'k')
xlim(ax1,[-1000 1000]);
ylim(ax1,[-1000 1000]);
xlim(ax2,[-1000 1000]);
ylim(ax2,[-1000 1000]);

xlabel(ax1,'true z (nm)')
ylabel(ax1,'fitted z (nm)')

xlabel(ax2,'true z (nm)')
ylabel(ax2,'corrected z (nm)')
end


function [Zint]=getZinterp(beads,Zintold,p)

% combine all coordiantes into large vectors to process them together
zobjectiveall=[];z0relativeall=[];zfitall=[];idall=[];dzall=[];
for k=1:length(beads)
    zobjectiveall=double(vertcat(zobjectiveall,beads(k).loc.zobjective));
    z0relativeall=double(vertcat(z0relativeall,beads(k).loc.zobjectiverelative)); %only for range
    zfitall=double(vertcat(zfitall,beads(k).loc.z));
    idall=double(vertcat(idall,k*ones(length(beads(k).loc.zobjective),1)));
    dzall=double(vertcat(dzall,beads(k).loc.dzcorr));          
end

% do interpolation in range where there are sufficient data points.
% qzfit=myquantile(zfitall,[0.05,0.95]);
% qzfit(1)=qzfit(1)+p.dz;qzfit(2)=qzfit(2)-p.dz;
inz=abs(z0relativeall)<p.maxrange; 
qzfit=[-1 1]*p.maxrange;
inz=inz&(zfitall)<qzfit(2)&(zfitall)>qzfit(1);
inz=inz&abs(dzall)<p.maxrange;
if ~isempty(Zintold)  %if next iteration: use last interpolation to identify and remove outliers
    dz=Zintold.interp(zobjectiveall,zfitall)-dzall; %distance from interpolation
    inz=inz&abs(dz)<p.cutoffrefine; % remove those data points
    h=histcounts(idall(inz),(1:max(idall)+1))';
    minpoints=p.maxrange/p.dz;
    innump=h(idall)>minpoints; %remove beads which do not have a minimum number of data points left
    inz=inz&innump;
end

%usually there are many beads close to the glass, this can lead to some
%quite rapid changes of the interpolation in the vicinity. To reduce this
%effect, we add additional anchor points for the interpolation to enzure
%dz==0 for small objective positions where no acchor points exist.
zfitx0=(qzfit(1):p.dz:qzfit(2))';
zfitx=repmat(zfitx0,10,1);
zfitallh=vertcat(zfitall(inz),zfitx,zfitx,zfitx,zfitx);
zobjectiveallh=vertcat(zobjectiveall(inz),min(zfitx)-zfitx,min(zfitx)/2-zfitx,min(zfitx)/4-zfitx,min(zobjectiveall(inz))-0*zfitx);   
dzallh=vertcat(dzall(inz),0*zfitx,0*zfitx,0*zfitx,0*zfitx); 

xrange=round(min(zobjectiveall(inz))/100)*100:100:max(zobjectiveall(inz));
yrange=round(qzfit(1)/10)*10:10: qzfit(2);
[X,Y]=meshgrid(xrange,yrange);  

%interpolation
Z=RegularizeData3D(zobjectiveallh,zfitallh,dzallh,xrange,yrange,'smoothness',[p.smoothframe p.smoothz],'extend','always');
Zint.interp=griddedInterpolant(X',Y',Z');

%plot result
if ~isempty(p.axhere)
    Zplot=Zint.interp(X',Y')';
    scatter3(p.axhere,zobjectiveall(inz),zfitall(inz),dzall(inz),1,'k')
    xlabel(p.axhere,'objective position above glass (nm)');ylabel(p.axhere,'zfit (nm)'); zlabel(p.axhere,'correction (nm)');
    hold(p.axhere,'on')
    Zcolor=Zplot; Zcolor(Zcolor>max(dzall(inz)))=max(dzall(inz));
    s=surf(p.axhere,X,Y,Zplot,Zcolor);
    s.FaceAlpha=0.8;
    s.EdgeColor='none';
    p.axhere.ZLim(1)=min(dzall(inz));
    p.axhere.ZLim(2)=max(dzall(inz));
    
end
end

function  f0glass=getf0glass(beads,p)
%determine the position of the glass as the robust minimum of bead
%positions
if isempty(p.axhere)
    f=figure;ax=gca;
else
    ax=p.axhere;
end

f0=[beads.f0];
dzh=50/p.dz;
induse=f0<dzh*60;
f0=f0(induse);    
range=min(f0):dzh:max(f0);
h=histogram(ax,f0,range);


[mh]=max(h.Values);
ind=find(h.Values>mh*.4,1,'first');
f0h=range(ind);
ind=(f0>f0h-2*dzh&f0<f0h+2*dzh);
f0glass=mean(f0(ind));
hold(ax, 'on');
xlabel('bead position (frame)')
ylabel('counts')
title('bead positions')

if isempty(p.axhere)
    close(f)
else
    plot(ax,f0glass,ones(size(f0glass)),'k*')
end
end

function err1=geterrors(beads,Zint)
yrange=Zint.interp.GridVectors{2};
for k=length(beads):-1:1
    zh=double(beads(k).loc.z);
    zglass=beads(k).loc.zobjective;
    z0f=beads(k).loc.dzcorr;
    inz=abs(zh<300) & abs(z0f)<300 & (zh)<yrange(end) & (zh)>yrange(1);
    dz=Zint.interp(zglass(inz),zh(inz))-z0f(inz);
    err1(k)=mean(dz.^2);
end
end

