function [spline,indgood,curves]=getspline_so(beads,p)
    curves=[];
    for B=length(beads):-1:1
        beadz0=(beads(B).f0)*p.dz;
        
        if contains(p.zcorr,'astig')||contains(p.zcorr,'corr')
             beadz=(beads(B).loc.frames*p.dz)-beadz0;
        else 
            beadz=(beads(B).loc.frames-p.midpoint)*p.dz;
        end
           sx=beads(B).loc.PSFxpix; 
           sy=beads(B).loc.PSFypix; 
           z=beadz; phot=beads(B).loc.phot;  
        inzr=z>=p.gaussrange(1)&z<=p.gaussrange(2);
        curves(B).sx=double(sx(inzr));
        curves(B).sy=double(sy(inzr));
        curves(B).z=double(z(inzr));
        curves(B).phot=double(phot(inzr));
        curves(B).xpos=beads(B).pos(1);
        curves(B).ypos=beads(B).pos(2);

    end
    
    %get calibrations
    bh=curves;
    p.ax=p.ax_z;
     [spline,indgood]=getcleanspline(curves,p);
     drawnow
end


function [s,indgood2]=getcleanspline(curves,p)

 za=vertcat(curves(:).z);
Sxa=vertcat(curves(:).sx);
Sya=vertcat(curves(:).sy);

indz=za>p.gaussrange(1)&za<p.gaussrange(2);
z=za(indz);
Sx=Sxa(indz);
Sy=Sya(indz);

warning('off','curvefit:fit:iterationLimitReached');
splinex=fit(z,Sx,'poly6','Robust','LAR','Normalize','on');
spliney=fit(z,Sy,'poly6','Robust','LAR','Normalize','on');
warning('on','curvefit:fit:iterationLimitReached');
hold off
for k=length(curves):-1:1
    w=(curves(k).phot);
    w=1;
    zh=curves(k).z;
    indzh=zh>p.gaussrange(1)&zh<p.gaussrange(2);
    w=w.*indzh;
    errh=(curves(k).sx-splinex(zh)).^2.*w+(curves(k).sy-spliney(zh)).^2.*w;
    err(k)=sqrt(sum(errh)/sum(w));
    errh2=(curves(k).sx-splinex(zh)).*w./curves(k).sx;
    errh3=(curves(k).sy-spliney(zh)).*w./curves(k).sy;
    err2(k)=abs(sum(errh2)/sum(w));
    err3(k)=abs(sum(errh3)/sum(w));
end

%  
[em,es]=robustMean(err);
if isnan(es), es=em; end
indgood2=err<em+2.5*es & err2+err3<.2;
if sum(indgood2)==0
    indgood2=true(1,length(curves));
end
zg=vertcat(curves(indgood2).z);
indz=zg>p.gaussrange(1)&zg<p.gaussrange(2);
zg=zg(indz);
sxg=vertcat(curves(indgood2).sx);
syg=vertcat(curves(indgood2).sy);
sxg=sxg(indz);
syg=syg(indz);

splinex2=getspline(sxg,zg,1./(abs(sxg-splinex(zg))+.1));
spliney2=getspline(syg,zg,1./(abs(syg-spliney(zg))+.1));

zt=min(zg):0.01:max(zg);
z1a=[];z2a=[];x1a=[];x2a=[];

for k=1:length(curves)
    if indgood2(k)
        z1a=[z1a; curves(k).z; curves(k).z];
        x1a=[x1a; curves(k).sx; curves(k).sy];
    else
        z2a=[z2a; curves(k).z; curves(k).z];
        x2a=[x2a; curves(k).sx; curves(k).sy];
    end
end
    if isempty(z2a) %only good curves: still plot a bad point in order to have legend correct.
        z2a=0;x2a=0;
    end
    plot(p.ax,z2a,x2a,'g.')
    hold(p.ax ,'on')
    plot(p.ax,z1a,x1a,'r.')
    plot(p.ax,zt,splinex2(zt),'k')
    plot(p.ax,zt,spliney2(zt),'k')

    xlim(p.ax,[zt(1) zt(end)])
    ylim(p.ax,[0 min(5,max(max(splinex2(zt)),max(spliney2(zt))))])
    xlabel(p.ax,'z (nm)');ylabel(p.ax,'PSFx, PSFy (pixel)');
    drawnow
    s.x=splinex2;s.y=spliney2;
    s.zrange=[zt(1) zt(end)];
    title(p.ax,'Lateral size of the PSF');


zr=s.zrange(1):1:s.zrange(2);
midp=round(length(zr)/8);


[~,ind1x]=max(s.x(zr(1:midp)));
[~,ind1y]=max(s.y(zr(1:midp)));
[~,ind2x]=max(s.x(zr(midp:end)));
[~,ind2y]=max(s.y(zr(midp:end)));

z1=max(zr(ind1x),zr(ind1y));
z2=min(zr(ind2x+midp-1),zr(ind2y+midp-1));

s.maxmaxrange=[z1 z2];
end

function spline=getspline(S,z,w,p)
if nargin<4
    p=0.96;
end
[zs,zind]=sort(z);Ss=S(zind);ws=w(zind);
spline=fit(zs,Ss,'smoothingspline','Weights',ws,'Normalize','on','SmoothingParam',p); 
end

