function transform=transform_locsN(locData,p)
%get localizations
if p.uselayers
    [locref,xx]=locData.getloc({'xnm','ynm','frame','filenumber','znm','PSFxnm'},'layer',p.reflayer.Value,'position','roi','grouping','ungrouped');
    [loctarget,yy]=locData.getloc({'xnm','ynm','frame','filenumber','znm','PSFxnm'},'layer',p.targetlayer.Value,'position','roi','grouping','ungrouped');  
else
    locref=locData.getloc({'xnm','ynm','frame','filenumber','inungrouped','znm','PSFxnm'},'position','roi');
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
    fn=fieldnames(locref);
    for k=1:length(fn)
        if ~isempty(locref.(fn{k}))
            locref.(fn{k})=locref.(fn{k})(ind);
        end
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

roinm=p.currentfileinfo.roi;
roinm([1 3])=roinm([1 3])*p.currentfileinfo.cam_pixelsize_um(1)*1000;
roinm([2 4])=roinm([2 4])*p.currentfileinfo.cam_pixelsize_um(end)*1000;

roiend(1)=p.currentfileinfo.roi(1)+p.currentfileinfo.roi(3);
roiend(2)=p.currentfileinfo.roi(2)+p.currentfileinfo.roi(4);
separator=roinm(1:2)+roinm(3:4);

if ~isfield(p,'repetition')
    p.repetition=[];
end
if ~p.useT %no initial transformation: use cross-correlation
    dx=0;dy=0;
    %set parameters for cross-correlation: where on the chip is reference,
    %target
switch p.targetpos.selection
    case 'top'
        separator(2)=roinm(2)+roinm(4)/2;
        spix=separator/(p.currentfileinfo.cam_pixelsize_um(1)*1000);
        indtarget=indtarget&loctarget.ynm<separator(2);
        indref=indref&locref.ynm>=separator(2);
        xrangecamr=[0 roiend(1)];yrangecamr=[spix(2) roiend(2)];
        xrangecamt=[0 roiend(1)];yrangecamt=[0 spix(2)]; 
        dy=separator(2);
    case 'bottom'
        separator(2)=roinm(2)+roinm(4)/2;
        spix=separator/(p.currentfileinfo.cam_pixelsize_um(1)*1000);
        indtarget=indtarget&loctarget.ynm>separator(2);
        indref=indref&locref.ynm<=separator(2);
        xrangecamr=[0 spix(1)];yrangecamr=[0 spix(2)];
        xrangecamt=[0 spix(1)];yrangecamt=[spix(2) roiend(2)];
        dy=-separator(2);
    case 'left'
        separator(1)=roinm(1)+roinm(3)/2;
        indtarget=indtarget&loctarget.xnm<separator(1);
        indref=indref&locref.xnm>=separator(1);
        spix=separator/(p.currentfileinfo.cam_pixelsize_um(1)*1000);
        xrangecamr=[spix(1) roiend(1)];yrangecamr=[0 spix(2)];
        xrangecamt=[0 spix(1)];yrangecamt=[0 spix(2)];
         dx=separator(1);
    case 'right'
        separator(1)=roinm(1)+roinm(3)/2;
        indtarget=loctarget.xnm>separator(1);
        indtarget=indtarget&loctarget.xnm>separator(1);
        indref=indref&locref.xnm<=separator(1);
        spix=separator/(p.currentfileinfo.cam_pixelsize_um(1)*1000);
        xrangecamr=[0 spix(1)];yrangecamr=[0 spix(2)];
        xrangecamt=[spix(1) roiend(1)];yrangecamt=[0 spix(2)];
        dx=-separator(1); %XXXXXX
        dy=0;
    otherwise
        dx=0;
        dy=0;
        fileref=locData.files.file(locref.filenumber(1));filetar=locData.files.file(loctarget.filenumber(1));
        xrangecamr=[0 fileref.info.roi(3)+fileref.info.roi(1)];
        yrangecamr=[0 fileref.info.roi(4)+fileref.info.roi(2)];
        xrangecamt=[0 filetar.info.roi(3)+filetar.info.roi(1)];
        yrangecamt=[0 filetar.info.roi(4)+filetar.info.roi(2)];
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
else %all initial estimation:
    loctT.x=loctT.xnm-0*dx+p.register_parameters.initialshiftx;
    loctT.y=loctT.ynm-0*dy+p.register_parameters.initialshifty;
    %take into account the mirroring
    midmirror=separator;
    switch p.targetmirror.selection
        case 'left-right'
            loctT.x=2*midmirror(1)-loctT.x;
            mirrorch2= 1;
        case 'up-down' 
            loctT.y=2*midmirror(2)-loctT.y;
            mirrorch2= 2;
        case 'both'  
            loctT.x=2*midmirror(1)-loctT.x;
            loctT.y=2*midmirror(2)-loctT.y;
            mirrorch2=[1 2];
        otherwise
            loctT.x=loctT.x+dx;
            loctT.y=loctT.y+dy;
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

if cutout
    loctT.frame(~indtarget)=-1;
end

pixelsizerec=p.register_parameters.pixelsizenm;
rangex=[min(min(locref.x),min(loctT.x)) max(max(locref.x),max(loctT.x))] +[-1 1]*pixelsizerec*5;
rangey=[min(min(locref.y),min(loctT.y)) max(max(locref.y),max(loctT.y))] +[-1 1]*pixelsizerec*5;

imr=myhist2(locref.x,locref.y,pixelsizerec,pixelsizerec,rangex,rangey);
imt=myhist2(loctT.x,loctT.y,pixelsizerec,pixelsizerec,rangex,rangey);

qimr=quantile(imr(:),.998);
if qimr>0
    imr(imr>qimr)=qimr;
end
qimt=quantile(imt(:),.998);
if qimt>0
    imt(imt>qimt)=qimt;
end

imr=sqrt(imr);
imt=sqrt(imt);

axim=initaxis(p.resultstabgroup,['img' p.repetition]);
imagesc(axim,horzcat(imr,imt))
title('images used for cross-correlation')
axis(axim, 'off')

maxshift=round(p.register_parameters.maxshift_corr/pixelsizerec);
ax1=initaxis(p.resultstabgroup,['corr' p.repetition]);

[dxpt,dypt]=getShiftCorr(imr,imt,1,maxshift);
dxcorr=dxpt*pixelsizerec;
dycorr=dypt*pixelsizerec;

loctT.x=loctT.x+dxcorr;
loctT.y=loctT.y+dycorr;

[iAa,iBa,na,nb,nseen]=matchlocsall(locref,loctT,0,0,p.register_parameters.maxshift_match,p.register_parameters.maxlocsused);


% x,y now in pixels!
% roitarget=locData.files.file(loctarget.filenumber(1)).info.roi;
% roiref=locData.files.file(locref.filenumber(1)).info.roi;
pixtarget=locData.files.file(loctarget.filenumber(1)).info.cam_pixelsize_um*1000;
pixref=locData.files.file(locref.filenumber(1)).info.cam_pixelsize_um*1000;


%still use camera chip as reference, not ROI
locref.x=locref.xnm/pixref(1);
locref.y=locref.ynm/pixref(end);
loctarget.x=loctarget.xnm/pixtarget(1);
loctarget.y=loctarget.ynm/pixtarget(end);

if p.useT
    transform=Tinitial.copy;
else
    transform=interfaces.LocTransformN;
    transform.setTransform(1,'type',p.transform.selection,'unit','pixel','parameter',p.transformparam,'cam_pixnm',pixref,'xrange',xrangecamr,'yrange',yrangecamr,'mirror',0);
    transform.setTransform(2,'type',p.transform.selection,'unit','pixel','parameter',p.transformparam,'cam_pixnm',pixtarget,'xrange',xrangecamt,'yrange',yrangecamt,'mirror',mirrorch2);
end

%XXXXXXX still need to include mirroring...XXXXXX
ztransform=isfield(locref,'znm')&&~isempty(locref.znm)&& any(locref.znm~=0) && any(loctarget.znm~=0) ;
if ztransform
    lref=horzcat(locref.x(iAa),locref.y(iAa),locref.znm(iAa));
    lt=horzcat(loctarget.x(iBa),loctarget.y(iBa),loctarget.znm(iBa));
else
    lref=horzcat(locref.x(iAa),locref.y(iAa));
    lt=horzcat(loctarget.x(iBa),loctarget.y(iBa));
end

transform.findTransform(2,lref,lt)

ltr=transform.transformToReference(2,lt);
   dx=ltr(:,1)-locref.x(iAa);
   dy=ltr(:,2)-locref.y(iAa);
   
   if ztransform
       dz=ltr(:,3)-locref.znm(iAa);
       dzb=loctarget.znm(iBa)-locref.znm(iAa);
       initaxis(p.resultstabgroup,['z' p.repetition])
       histogram(dz);hold on ; histogram(dzb);hold off
       xlabel('z_{ref}-z_{target} (nm)')
   end
   
initaxis(p.resultstabgroup,['dxy' p.repetition])
   dscatter(dx,dy)
   xlabel('x_{ref}-x_{target} (camera pixels)')
   ylabel('y_{ref}-y_{target} (camera pixels)')

   title({['number of anchor points: ' num2str(length(iBa)) ' of ' num2str(nseen) '=' num2str(length(iBa)/nseen*100,'%2.1f') '%']...
       ,['dx= ' num2str(std(dx),3) ' pix, dy= ' num2str(std(dy),3) ' pix']});

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
end
