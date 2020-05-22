function transform=transform_locsN(locData,p)

%get fields
if p.uselayers
    [locref,xx]=locData.getloc({'xnm','ynm','frame','filenumber','znm','PSFxnm'},'layer',p.reflayer.Value,'position','roi','grouping','ungrouped');
    [loctarget,yy]=locData.getloc({'xnm','ynm','frame','filenumber','znm','PSFxnm'},'layer',p.targetlayer.Value,'position','roi','grouping','ungrouped');
    
else
    locref=locData.getloc({'xnm','ynm','frame','filenumber','inungrouped','znm','PSFxnm'},'position','roi');
%     loctarget=locData.getloc('xnm','ynm','frame','filenumber','znm','PSFxnm');
    
    filter=locData.layer(1).filter;
    ind=true(length(locref.xnm),1);
    if isfield(filter,'locprecnm')
        ind=ind&filter.locprecnm(locref.inungrouped);
    end
    if isfield(filter,'PSFxnm')
        ind=ind&filter.PSFxnm(locref.inungrouped);
    end
    if isfield(filter,'frame')
        ind=ind&filter.frame(locref.inungrouped);
    end
    if ~p.allfiles
    ind=ind&locref.filenumber==p.dataselect.Value;
    end
%     indref=ind;
%     indtarget=ind;
    fn=fieldnames(locref);
    for k=1:length(fn)
        if ~isempty(locref.(fn{k}))
        locref.(fn{k})=locref.(fn{k})(ind);
        end
%         loctarget.(fn{k})=loctarget.(fn{k})(ind);
    end
    loctarget=locref;
end
indtarget=true(size(loctarget.frame));
indref=true(size(locref.frame));
locref.x=locref.xnm;
locref.y=locref.ynm;
loctarget.x=loctarget.xnm;
loctarget.y=loctarget.ynm;
loctT=loctarget;
% file=locData.files.file(p.layer1_.ch_filelist.Value);

roinm=p.currentfileinfo.roi;
roinm([1 3])=roinm([1 3])*p.currentfileinfo.cam_pixelsize_um(1)*1000;
roinm([2 4])=roinm([2 4])*p.currentfileinfo.cam_pixelsize_um(end)*1000;

roiend(1)=p.currentfileinfo.roi(1)+p.currentfileinfo.roi(3);
roiend(2)=p.currentfileinfo.roi(2)+p.currentfileinfo.roi(4);
chipsizenm=p.currentfileinfo.cam_pixelsize_um*1000.*[p.currentfileinfo.Width p.currentfileinfo.Height]; 
% facsize=ones(2,1);
% separator=chipsizenm;
separator=roinm(1:2)+roinm(3:4);
spix=separator/(p.currentfileinfo.cam_pixelsize_um(1)*1000)/2;
mirroradd=zeros(2,1);

if ~isfield(p,'repetition')
    p.repetition=[];
end
if p.useT
else
switch p.targetpos.selection
    case 'top'
%         dy=-chipsizenm(2)/2;
        dy=0;
        dx=0;
%         separator(2)=chipsizenm(2)/2;
        separator(2)=roinm(2)+roinm(4)/2;
        spix=separator/(p.currentfileinfo.cam_pixelsize_um(1)*1000);
        indtarget=indtarget&loctarget.ynm<separator(2);
        indref=indref&locref.ynm>=separator(2);
        xrangecamr=[0 roiend(1)];yrangecamr=[spix(2) roiend(2)];
        xrangecamt=[0 roiend(1)];yrangecamt=[0 spix(2)];
        mirroradd(2)=chipsizenm(2);
    case 'bottom'
%         dy=chipsizenm(2)/2;
        dy=0;
        dx=0;
%         separator(2)=chipsizenm(2)/2;
        separator(2)=roinm(2)+roinm(4)/2;
        spix=separator/(p.currentfileinfo.cam_pixelsize_um(1)*1000);
        indtarget=indtarget&loctarget.ynm>separator(2);
        indref=indref&locref.ynm<=separator(2);
        xrangecamr=[0 spix(1)];yrangecamr=[0 spix(2)];
        xrangecamt=[0 spix(1)];yrangecamt=[spix(2) roiend(2)];
    case 'left'
%         dx=-chipsizenm(1)/2;
        dx=0;
        dy=0;
%         separator(1)=chipsizenm(1)/2;
        separator(1)=roinm(1)+roinm(3)/2;
        indtarget=indtarget&loctarget.xnm<separator(1);
        indref=indref&locref.xnm>=separator(1);
        spix=separator/(p.currentfileinfo.cam_pixelsize_um(1)*1000);
        xrangecamr=[spix(1) roiend(1)];yrangecamr=[0 spix(2)];
        xrangecamt=[0 spix(1)];yrangecamt=[0 spix(2)];
         mirroradd(1)=chipsizenm(2);
    case 'right'
%         dx=chipsizenm(1)/2;
        dx=0; %XXXXXX
        dy=0;
%         separator(1)=chipsizenm(1)/2;
        separator(1)=roinm(1)+roinm(3)/2;
        indtarget=loctarget.xnm>separator(1);
        indtarget=indtarget&loctarget.xnm>separator(1);
        indref=indref&locref.xnm<=separator(1);
        spix=separator/(p.currentfileinfo.cam_pixelsize_um(1)*1000);
        xrangecamr=[0 spix(1)];yrangecamr=[0 spix(2)];
        xrangecamt=[spix(1) roiend(1)];yrangecamt=[0 spix(2)];
    otherwise
        dx=0;
        dy=0;
        fileref=locData.files.file(locref.filenumber(1));filetar=locData.files.file(loctarget.filenumber(1));
        xrangecamr=[0 fileref.info.roi(3)+fileref.info.roi(1)];
        yrangecamr=[0 fileref.info.roi(4)+fileref.info.roi(2)];
        xrangecamt=[0 filetar.info.roi(3)+filetar.info.roi(1)];
        yrangecamt=[0 filetar.info.roi(4)+filetar.info.roi(2)];
%         indtarget=true(size(loctarget.xnm));
end
end

cutout=false;
if p.useT
    if ~ischar(p.Tfile)
        Tinitial=p.Tfile;
    else
        Tinitial=loadtransformation(locData,p.Tfile,p.dataselect.Value);
    end
    if isa(Tinitial,'interfaces.LocTransformN')
%         Tinitial.setTransform([],'cam_pixnm',p.currentfileinfo.cam_pixelsize_um*1000)
        pos=Tinitial.transformToReference(2,horzcat(loctarget.xnm,loctarget.ynm),'nm');
        loctT.x=pos(:,1);loctT.y=pos(:,2);
           indref=Tinitial.getPart(1,horzcat(locref.xnm,locref.ynm),'nm');
           indtarget=Tinitial.getPart(2,horzcat(loctarget.xnm,loctarget.ynm),'nm');
    else
        [loctT.x,loctT.y]=Tinitial.transformCoordinatesInv(loctarget.xnm,loctarget.ynm);
        mirrorinfo=Tinitial.tinfo.mirror;
        if contains(mirrorinfo.targetmirror,'no')
            cutout=true;
        end
    end
    %     pos=Tinitial.pos;
%     size=Tinitial.size;
else %all initial estimation:
%     approximate shift from size and position
   
    loctT.x=loctT.xnm-dx+p.register_parameters.initialshiftx;
    loctT.y=loctT.ynm-dy+p.register_parameters.initialshifty;
    
    midmirror=separator;
    %  initial mirror
%     qx=myquantile(vertcat(locref.x,loctarget.x),[0.01,0.99]);
%     midmirror(1)=mean(qx);
%     qy=myquantile(vertcat(locref.y,loctarget.y),[0.01,0.99]);
%     midmirror(2)=mean(qy);
%     midmirror=chipsizenm/4;
%     mirrorinfo.midmirror=midmirror;
%     mirrorinfo.targetmirror=p.targetmirror.selection;
    
    switch p.targetmirror.selection
        case 'left-right'
            loctT.x=2*midmirror(1)-loctT.x+mirroradd(1);
            mirrorch2= 1;
        case 'up-down' 
            loctT.y=2*midmirror(2)-loctT.y+mirroradd(2);
            mirrorch2= 2;
        case 'both'  
            loctT.x=2*midmirror(1)-loctT.x+mirroradd(1);
            loctT.y=2*midmirror(2)-loctT.y+mirroradd(2);
            mirrorch2=[1 2];
        otherwise
            cutout=true;
            mirrorch2=0;
           
    end
    if isfield(p.register_parameters,'initial_mag') && p.register_parameters.initial_mag~=1
        loctT.x=loctT.x*p.register_parameters.initial_mag;
        loctT.y=loctT.y*p.register_parameters.initial_mag;
    end

end

%reduce
locref=copystructReduce(locref,indref);
loctarget=copystructReduce(loctarget,indtarget);
loctT=copystructReduce(loctT,indtarget);

% indref=loctarget.xnm<chipsizenm*facsize(1)&loctarget.ynm<chipsizenm*facsize(2);
% loctT=copystructReduce(loctT,ind);
if cutout
loctT.frame(~indtarget)=-1;
end

pixelsizerec=p.register_parameters.pixelsizenm;
% if p.uselayers
    rangex=[min(min(locref.x),min(loctT.x)) max(max(locref.x),max(loctT.x))] +[-1 1]*pixelsizerec*5;
    rangey=[min(min(locref.y),min(loctT.y)) max(max(locref.y),max(loctT.y))] +[-1 1]*pixelsizerec*5;
% else
%     
% roi=p.currentfileinfo.roi;
% roinm=roi;
% roinm([1 3])=roinm([1 3])*p.currentfileinfo.cam_pixelsize_um(1)*1000;
% roinm([2 4])=roinm([2 4])*p.currentfileinfo.cam_pixelsize_um(end)*1000;
% 
% 
% % pos=[mean(roinm([1,3])) meannm(roi([2,4]))];
% %     rsize=[roi(3)-roi(1) roi(4)-roi(2)];
% 
% rangex=[roinm(1) roinm(1)+roinm(3)*facsize(1)];
% rangey=[roinm(2) roinm(2)+roinm(4)*facsize(2)];
% end

% imr=myhist2(locref.x(indref),locref.y(indref),pixelsizerec,pixelsizerec,rangex,rangey);
% imt=myhist2(loctT.x(indtarget),loctT.y(indtarget),pixelsizerec,pixelsizerec,rangex,rangey);

imr=myhist2(locref.x,locref.y,pixelsizerec,pixelsizerec,rangex,rangey);
imt=myhist2(loctT.x,loctT.y,pixelsizerec,pixelsizerec,rangex,rangey);

qimr=quantile(imr(:),.998);
imr(imr>qimr)=qimr;
qimt=quantile(imt(:),.998);
imt(imt>qimt)=qimt;

imr=sqrt(imr);
imt=sqrt(imt);

axim=initaxis(p.resultstabgroup,['img' p.repetition]);
imagesc(axim,horzcat(imr,imt))
title('images used for cross-correlation')
axis(axim, 'off')
% implot=zeros(size(imr,1),size(imr,2),3);
% implot(:,:,1)=imr;implot(:,:,2)=imt;
% ax1=initaxis(p.resultstabgroup,'images');
% 
% figure(99);
% hold off
% plot(locref.x,locref.y,'bo')
% hold on
% plot(loctT.x,loctT.y,'rx')
% % plot(loctarget.xnm,loctarget.ynm,'g.')
% 
% 


maxshift=round(p.register_parameters.maxshift_corr/pixelsizerec);
ax1=initaxis(p.resultstabgroup,['corr' p.repetition]);

% s=size(imr)
% ima=zeros(s(1),s(2),3);
% ima(:,:,1)=imr;ima(:,:,2)=imt;
% figure(88);imagesc(ima)
% asdf
[dxpt,dypt]=getShiftCorr(imr,imt,1,maxshift);
dxcorr=dxpt*pixelsizerec;
dycorr=dypt*pixelsizerec;

loctT.x=loctT.x+dxcorr;
loctT.y=loctT.y+dycorr;
% [loctarget,indt]=transform.getPartLocs(loc,loc.x(ig),loc.y(ig),'target');
% [locref,indr]=transform.getPartLocs(loc,loc.x(ig),loc.y(ig),'reference');
% figure(99)
% plot(loctT.x(indtarget),loctT.y(indtarget),'gx',locref.x,locref.y,'bo')
% hold off

% locrefred=reduceFieldLength(locref,ig);
% loctargetred=reduceFieldLength(loctarget,ig);

[iAa,iBa,na,nb,nseen]=matchlocsall(locref,loctT,0,0,p.register_parameters.maxshift_match,p.register_parameters.maxlocsused);


% x,y now in pixels!
roitarget=locData.files.file(loctarget.filenumber(1)).info.roi;
roiref=locData.files.file(locref.filenumber(1)).info.roi;
pixtarget=locData.files.file(loctarget.filenumber(1)).info.cam_pixelsize_um*1000;
pixref=locData.files.file(locref.filenumber(1)).info.cam_pixelsize_um*1000;

% locref.x=locref.xnm/pixref(1)-roiref(1);
% locref.y=locref.ynm/pixref(2)-roiref(2);
% loctarget.x=loctarget.xnm/pixtarget(1)-roitarget(1);
% loctarget.y=loctarget.ynm/pixtarget(2)-roitarget(2);

%still use camera chip as reference, not ROI
locref.x=locref.xnm/pixref(1);
locref.y=locref.ynm/pixref(end);
loctarget.x=loctarget.xnm/pixtarget(1);
loctarget.y=loctarget.ynm/pixtarget(end);

% transform.findTransform(locref.x(iAa),locref.y(iAa),loctarget.x(iBa),loctarget.y(iBa))

if p.useT
    transform=Tinitial.copy;
else
%     if contains(mirrorinfo.targetmirror,'up-down')
%         mirrorch2= 2;
%     elseif contains(mirrorinfo.targetmirror,'left-right')
%         mirrorch2= 1;
%     else
%         mirrorch2=0;
%     end
    transform=interfaces.LocTransformN;
    transform.setTransform(1,'type',p.transform.selection,'unit','pixel','parameter',p.transformparam,'cam_pixnm',pixref,'xrange',xrangecamr,'yrange',yrangecamr,'mirror',0);
    transform.setTransform(2,'type',p.transform.selection,'unit','pixel','parameter',p.transformparam,'cam_pixnm',pixtarget,'xrange',xrangecamt,'yrange',yrangecamt,'mirror',mirrorch2);
end

%XXXXXXX still need to include mirroring...XXXXXX
% t.parameter=p.transformparam;




ztransform=isfield(locref,'znm')&&~isempty(locref.znm);
if ztransform
    lref=horzcat(locref.x(iAa),locref.y(iAa),locref.znm(iAa));
    lt=horzcat(loctarget.x(iBa),loctarget.y(iBa),loctarget.znm(iBa));
else
    lref=horzcat(locref.x(iAa),locref.y(iAa));
    lt=horzcat(loctarget.x(iBa),loctarget.y(iBa));
end

transform.findTransform(2,lref,lt)

ltr=transform.transformToReference(2,lt);
% transform.findTransform(locref.x(iAa),locref.y(iAa),loctT.x(iBa),loctT.y(iBa),t)
% if p.showresults
    
   
%         [xa, ya]=transform.transformCoordinatesInv((loctarget.x(iBa)),(loctarget.y(iBa)));
%      [xa, ya]=transform.transformCoordinatesInv((loctT.x(iBa)),(loctT.y(iBa)));
%  [xa, ya]=transform.transformCoordinates((loc.x(indr(iBa))),(loc.y(indt(iBa))),'target');
 
   dx=ltr(:,1)-locref.x(iAa);
   dy=ltr(:,2)-locref.y(iAa);
   
   if ztransform
       dz=ltr(:,3)-locref.znm(iAa);
       dzb=loctarget.znm(iBa)-locref.znm(iAa);
       initaxis(p.resultstabgroup,['z' p.repetition])
       histogram(dz);hold on ; histogram(dzb);hold off
       xlabel('z_{ref}-z_{target} (nm)')
   end
   
%    dx=loctarget.x(iBa)-locref.x(iAa);
%    dx=loctarget.y(iBa)-locref.y(iAa);
%    figure(88)
initaxis(p.resultstabgroup,['dxy' p.repetition])
   dscatter(dx,dy)
   xlabel('x_{ref}-x_{target} (camera pixels)')
   ylabel('y_{ref}-y_{target} (camera pixels)')

   title({['number of anchor points: ' num2str(length(iBa)) ' of ' num2str(nseen) '=' num2str(length(iBa)/nseen*100,'%2.1f') '%']...
       ,['dx= ' num2str(std(dx),3) ' pix, dy= ' num2str(std(dy),3) ' pix']});
%  ax3=initaxis(p.resultstabgroup,['hist' p.repetition]);
%  hist(dx,50)
%  xlabel('x_{ref}-x_{target} (camera pixels)')
 rr=rand(2000,1);
 ra=ceil(rr*length(iBa));
  rb=ceil(rr*length(nb));
  if isempty(nb)
      rb=[];
  end
 ax4=initaxis(p.resultstabgroup,['pos' p.repetition]);
 plot(loctarget.x(nb(rb)),loctarget.y(nb(rb)),'r+',loctarget.x(iBa(ra)),loctarget.y(iBa(ra)),'bx')
 legend('unpaired','paired');
 title({'random subset of localizations.';'Paired ones should span entire field'})
    xlabel('x (camera pixels)')
   ylabel('y (camera pixels)')
% transform.tinfo.targetpos=p.targetpos.selection;
% transform.tinfo.separator=separator;
% transform.tinfo.mirror=mirrorinfo;
% transform.tinfo.cam_pixelsize_nm=p.currentfileinfo.cam_pixelsize_um*1000;

