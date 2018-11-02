function drift=driftcorrection3D_so(x,y,z,frame,p)
pard.correctxy=true; %correct in x, y
pard.drift_timepoints=10; %number of time points (blocks) in x,y,
pard.drift_pixrec=5; %size of pixels for reconstruction (in units of x,y, e.g. nm)
pard.drift_window=7; %Size of window in which maximum of cross-correlation is found (in reconstructed pixels).
pard.drift_maxdrift=1000; %maximum drift (in units of x, eg.g nm)
pard.drift_maxpixels=4096; %maximum size of reconstructed image (to avoid memory overflow).

pard.correctz=true; % correct in z;
pard.drift_timepointsz=20; %time points in z
pard.drift_pixrecz=5; %pixelsize for reconstruction in z (nm)
pard.drift_windowz=9; %size of window used for maximum finding
pard.zrange=[-300 300]; %range in z used for drift correction;
pard.slicewidth=200;% width of slice for z-correlation (in units of x, eg. nm)

pard.smoothmode.Value=2; %one of: '1: smoothing cubic spline','2: linear'
pard.smoothpar=[]; %smoothing parameter for cubic spline. Empty for automatic.
pard.drift_reference=false; %set true if reference is last frame;
pard.drift_mirror2c.Value=1; %one of: 1: 'no mirror', 2: '2 Channels, mirrored, vertical split', 3: '2 Channels, mirrored, horizontal split'
pard.showresults=true;

pard.maxframeall=max(frame); %range in frames (time) used for correction. Outside uses edge values
pard.framestart=min(frame);
pard.framestop=pard.maxframeall;
p=copyfields(pard,p);

p.roi=[min(x),min(y),max(x)-min(x),max(y)-min(y)]; %roi used (x,y,wx,wy)

f=figure;
p.resultstabgroup=uitabgroup(f);
if ~p.correctxy && ~p.correctz
    return
end
locs.xnm=x;locs.ynm=y;locs.znm=z;locs.frame=frame;
[drift,driftinfo,fieldc]=getxyzdrift(locs,p);        
% locsnew=applydriftcorrection(drift,locs);
end


function [drift,driftinfo,fieldc]=getxyzdrift(locs,p)

drift=[];driftinfo=[];
                
if p.correctxy
   
    switch p.drift_mirror2c.Value
        case 1 %all
            ind=true(length(locs.xnm),1);
            rep=false;
            mirror='none';
            midpoint=0;
             p.repetitionname='';
        case 2 %horizontal
            midpoint=p.cam_pixelsize_nm(1)*(p.roi(1)+(p.roi(1)+p.roi(3))/2);
            ind=locs.xnm<=midpoint;
            rep=true;
            mirror='horizontal';
             p.repetitionname='1';
        case 3 %vertical
            midpoint=p.cam_pixelsize_nm(2)*(p.roi(2)+(p.roi(2)+p.roi(4))/2);
            ind=locs.ynm<=midpoint;
            rep=true;
            mirror='vertical';
            p.repetitionname='1';
    end
    [driftxy,driftinfoxy]=finddriftfeature(copystructReduce(locs,ind),p);
    driftinfoxy.mirror=mirror; driftinfoxy.midpoint=midpoint;
    driftinfo.xy=driftinfoxy;
    drift.xy=copyfields([],driftxy,{'x','y'});
    drift.xy(1).mirror=mirror;drift(1).xy(1).midpoint=midpoint;

    if rep
         p.repetitionname='2';
        [driftxy,driftinfoxy]=finddriftfeature(copystructReduce(locs,~ind),p);
        driftinfoxy.mirror=mirror; driftinfoxy.midpoint=midpoint;
        driftinfo.xy(2)=(driftinfoxy);
        drift.xy(2)=copyfields(drift.xy(1),driftxy,{'x','y'});
%         drift.xy(2).mirror=mirror;drift(2).midpoint=midpoint;
    end
end
if ~isempty(locs.znm)&&p.correctz
    if p.correctxy
        locsnew=copyfields(locs,applydriftcorrection(drift,locs),{'xnm','ynm'});
    else
        locsnew=locs;
        drift=[];
    end
    [driftz,driftinfoz]=finddriftfeatureZ(locsnew,p);
    
    drift.z=driftz.z;%copyfields(drift,driftz,'z');
     driftinfo.z=driftinfoz;%copyfields(driftinfo,driftinfoz);
     fieldc={'xnm','ynm','znm'};
else
    fieldc={'xnm','ynm'};
end
% locsall=copyfields([],obj.locData.loc,{fieldc{:},'frame','filenumber'});
% locsnew=applydriftcorrection(drift,locsall);
end

