function calibrate_globalworkflow(p)

if ~p.isglobalfit %only normal local calibratin, just call old proram
    calibrate3D_g(p);
    return
end

global S beadpos1 beadpos2
f=figure('Name','Bead calibration');
tg=uitabgroup(f);


if p.makeT || isempty(p.Tfile)
pr=getranges(p) ;
ph=p;
ph.isglobalfit=false;
%set spatial calibration
ph.outputfile={};
if contains(p.filelist{1},';') %each channel in separate file
    for k=1:length(p.filelist)
        ind=strfind(p.filelist{k},';');
        filelist1{k}=p.filelist{k}(1:ind-1);
        filelist2{k}=p.filelist{k}(ind+1:end);
    end
else
    filelist1=p.filelist; filelist2=p.filelist;
end


t1=uitab(tg,'Title','first channel');
ph.filelist=filelist1;
ph.tabgroup=  uitabgroup(t1);  
ph.yrange=pr.yrange1;ph.xrange=pr.xrange1;
ph.filechannel=1;
[S1,beadpos1,parameters1]=calibrate3D_g(ph);



t2=uitab(tg,'Title','second channel');
ph.filelist=filelist2;
ph.tabgroup=  uitabgroup(t2);  
ph.yrange=pr.yrange2;ph.xrange=pr.xrange2;
ph.filechannel=2;
[S2,beadpos2,parameters2]=calibrate3D_g(ph);


% [S,beadpos]=calibrate3D_g(ph);

% Later: also do test-fitting with corresponding spline coefficients
tt=uitab(tg,'Title','transformation');
p.tabgroup=  uitabgroup(tt);
% p.separator=p.Tsplitpos+parameters1.roi{1}(pr.roiind);
p.separator=p.Tsplitpos;
% find transform
% if p.makeT || isempty(p.Tfile)
%     transform=transform_locs_simple(beadpos1{1},beadpos2{1},p);
    transform=makeglobalTransform(beadpos1{1},beadpos2{1},p);
else
    l=load(p.Tfile);
    transform=l.transformation;
    p.Tmode=transform.tinfo.mirror.targetmirror;
    p.Tsplitpos=transform.tinfo.separator;
    pr=getranges(p) ;
%     split= transform.mirror;
end


ph=p;
ph.outputfile=[];
ph.isglobalfit=true;
ph.Tfile=transform;
% ph.outputfile=p.outputfile;
ph.outputfile=[];
t4=uitab(tg,'Title','global cal');
ph.tabgroup=  uitabgroup(t4);  
% ph.yrange=yrange1;ph.xrange=xrange1;
ph.filelist=p.filelist;
[S,beadpos,parameters_g]=calibrate3D_g(ph);

% if ~exist('S1','var') %take global one apart...
% recover S1, S2 dfrom S_glob to ensure same z reference, and PSF from
% corresponding beads.
    S1=S;%XXXXXX also take the right coefficients!!!
    S1.PSF=S1.PSF(1);
    S1.cspline.coeff={S1.cspline.global.coeffrawref};
    S1.cspline.normf=S1.cspline.normf(1);
    S1.cspline.mirror=0;
    S2=S;
    S2.PSF=S2.PSF(2);
    S2.cspline.coeff={S2.cspline.global.coeffrawtar};
    S2.cspline.normf=S2.cspline.normf(2);
    S1.Xrange=pr.xrange1;S2.Xrange=pr.xrange2;
    S1.Yrange=pr.yrange1;S2.Yrange=pr.yrange2;
% end
S1.Yrangeall=pr.yrangeall;S1.Xrangeall=pr.xrangeall;
S2.Yrangeall=pr.yrangeall;S2.Xrangeall=pr.xrangeall;
S2.posind=pr.XYpos;

if strcmp(pr.split,'rl') 
    SXY(1:length(S1),1)=S1;
    SXY(end+1:end+length(S2),1)=S2;
else
    SXY(1,1:length(S1))=S1;
    SXY(1,end+1:end+length(S2))=S2;
end
SXY_g=S;
transformation=parameters_g.transformation;
calibrationfigure=f;
if ~isempty(p.outputfile)
    if p.smap
        parameters1.smappos.P=[]; parameters2.smappos.P=[]; parameters_g.smappos.P=[];
        save(p.outputfile,'SXY','SXY_g','parameters_g','parameters1','parameters2','transformation');
        
    else
        save(p.outputfile,'gausscal','cspline_all','gauss_sx2_sy2','gauss_zfit','cspline','parameters');
    end
    filefig=strrep(p.outputfile,'.mat','.fig');
    savefig(calibrationfigure,filefig,'compact');
end
end



function pr=getranges(p)

if ~isfield(p,'yrange')
    p.yrange=[-inf inf];
end
if  ~isfield(p,'xrange')
    p.xrange=[-inf inf];
end

pr.mirror=contains(p.Tmode,'mirror');
switch p.Tmode
    case {'up-down','up-down mirror'}
        splitpos=p.Tsplitpos(1);
        if max(p.yrange)<splitpos %defined only in upper part
            yrange1=p.yrange;
            if contains(p.Tmode,'mirror')
                yrange2=sort(-p.yrange+2*splitpos);yrange2(yrange2<splitpos)=splitpos;
            else
                yrange2=sort(p.yrange+splitpos);yrange2(yrange2<splitpos)=splitpos;
            end
        else
            yrange1=([p.yrange splitpos]);yrange1(yrange1>splitpos)=splitpos;yrange1=unique(yrange1);
            yrange2=([p.yrange+ splitpos]);yrange2(yrange2<splitpos)=splitpos;yrange2=unique(yrange2);
        end
            
%         yrange=unique([p.yrange splitpos p.yrange+splitpos]);
%         yrange1=yrange(yrange<=splitpos);yrange2=yrange(yrange>=splitpos);
        pr.xrange1=p.xrange;pr.xrange2=p.xrange;
        yrangeall=[yrange1(1:end-1) splitpos yrange2(2:end)];
        pr.yrangeall=yrangeall;
        yrange1(end)=yrange1(end)-p.mindistance; %do not take into account locs too close to separator
        yrange2(1)=yrange2(1)+p.mindistance;
         pr.yrange1=yrange1;pr.yrange2=yrange2;
        pr.xrangeall=p.xrange;
        pr.XYpos=[1,2];
        
        pr.split='ud';
        pr.roiind=2;
    case {'right-left','right-left mirror'}
          splitpos=p.Tsplitpos(end);
         if max(p.xrange)<splitpos %defined only in upper part
            xrange1=p.xrange;
            if contains(p.Tmode,'mirror')
                xrange2=sort(-p.xrange+2*splitpos);xrange2(xrange2<splitpos)=splitpos;
            else
                xrange2=sort(p.xrange+splitpos);xrange2(xrange2<splitpos)=splitpos;
            end
        else
            xrange1=([p.xrange splitpos]);xrange1(xrange1>splitpos)=splitpos;xrange1=unique(xrange1);
            xrange2=([p.xrange+ splitpos]);xrange2(xrange2<splitpos)=splitpos;xrange2=unique(xrange2);
         end
            
        pr.yrange1=p.yrange;pr.yrange2=p.yrange;
        xrangeall=[xrange1(1:end-1) splitpos xrange2(2:end)];
        pr.xrangeall=xrangeall;
        
        xrange1(end)=xrange1(end)-p.mindistance; %do not take into account locs too close to separator
        xrange2(1)=xrange2(1)+p.mindistance;
        pr.xrange1=xrange1;
        pr.xrange2=xrange2;
        pr.yrangeall=p.yrange;
        pr.XYpos=[2,1];
        
        pr.split='rl';
        pr.roiind=1;
    case {'2 cam','2 cam u-d mirror','2 cam r-l mirror'}
        pr.xrange1=p.xrange;pr.yrange1=p.yrange;
        pr.xrange2=p.xrange;pr.yrange2=p.yrange;
        pr.roiind=1;
        pr.xrangeall=p.xrange; pr.yrangeall=p.yrange;
        pr.XYpos=[1,1];
        pr.split ='none';
end
%test: exchange channels
if p.switchchannels
xr1temp=pr.xrange1; 
pr.xrange1=pr.xrange2;
pr.xrange2=xr1temp;
yr1temp=pr.yrange1; 
pr.yrange1=pr.yrange2;
pr.yrange2=yr1temp;
end

end

% 
function transform=makeglobalTransform(bead1,bead2,ph)
%calculate transformN
pp=getranges(ph);
transform=interfaces.LocTransformN;
pt.mirror=[false false]; %ref
pt.xrange=pp.xrange1;
pt.yrange=pp.yrange1;
pt.unit='pixel';
pt.type='projective';
transform.setTransform(1,pt)
pt.mirror=0;

if pp.mirror
if contains(pp.split,'rl')
    pt.mirror= 1;
elseif contains(pp.split,'ud')
    pt.mirror= 2;
else
    pt.mirror=0;
end
end
pt.xrange=pp.xrange2;
pt.yrange=pp.yrange2;
transform.setTransform(2,pt)
th=ph.tabgroup;
ph.ax=th;
numf=5;
mp=round(size(bead1.x,1)+1)/2;
range=mp-numf:mp+numf;
bc1=horzcat(reshape(bead1.x(range,:),[],1),...
    reshape(bead1.y(range,:),[],1),...
    reshape(bead1.z(range,:),[],1),...
    reshape(bead1.filenumber(range,:)*100+bead1.frame(range,:),[],1));
bc2=horzcat(reshape(bead2.x(range,:),[],1),...
    reshape(bead2.y(range,:),[],1),...
    reshape(bead2.z(range,:),[],1),...
    reshape(bead2.filenumber(range,:)*100+bead2.frame(range,:),[],1));
ph.sepscale=5;
[transform ,iAa,iBa]=transform_locs_simpleN(transform,1, bc1,2,bc2,ph); 

end
% 

