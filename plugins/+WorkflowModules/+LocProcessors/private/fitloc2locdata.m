function locdat=fitloc2locdata(obj,locs,indin)
fieldsremove={'xerrpix','yerrpix','PSFxpix','PSFypix','xpix','ypix'};
fieldsremove={'xerrpix','yerrpix'};
keeptype={'frame','filenumber','inframes'};

if isfield(locs,'xerrpix') 
    locs.xpixerr=locs.xerrpix;
end
if isfield(locs,'yerrpix') 
    locs.ypixerr=locs.yerrpix;
end
fn=fieldnames(locs);
keepfields=setdiff(fn,fieldsremove);
for k=1:length(keepfields)
    if any(strcmp(keepfields{k},keeptype))
    locdat.(keepfields{k})=locs.(keepfields{k})(indin);
    else
        locdat.(keepfields{k})=single(locs.(keepfields{k})(indin));
    end
end
if isprop(obj,'fileinfo')
    pixelsize=obj.fileinfo.cam_pixelsize_um*1000;
    if myisfield(obj.fileinfo,'roi')&&~isempty(obj.fileinfo.roi)
    roi=obj.fileinfo.roi;
    else
        roi=zeros(2);
    end
else
    pixelsize=[1 1];
    roi=zeros(2);
end

locdat.xnm=(locs.xpix(indin)+roi(1))*pixelsize(1);
locdat.ynm=(locs.ypix(indin)+roi(2))*pixelsize(end);

locdat.xerrpix=locdat.xnm*0+1;
% if isfield(locs,'xerrpix')
% locdat.xerr=locs.xerrpix(indin)*pixelsize(1);
% end
% if isfield(locs,'xpixerr')
% locdat.xerr=locs.xpixerr(indin)*pixelsize(1);
% end
% if isfield(locs,'xpix2err')
% locdat.xerr2=locs.xpix2err(indin)*pixelsize(1);
% locdat.xerr1=locdat.xerr;
% locdat.xerr=min(locdat.xerr1,locdat.xerr2);
% end

locdat.yerrpix=locdat.xerrpix;
% if isfield(locs,'yerrpix')
% locdat.yerr=locs.yerrpix(indin)*pixelsize(end);
% end
% if isfield(locs,'ypixerr')
% locdat.yerr=locs.ypixerr(indin)*pixelsize(end);
% end
% if isfield(locs,'ypix2err')
% locdat.yerr2=locs.ypix2err(indin)*pixelsize(end);
% locdat.yerr1=locdat.yerr;
% locdat.yerr=min(locdat.yerr1,locdat.yerr2);
% end

% if isfield(locs,'xpix2')
% locdat.xnm2=(locs.xpix2(indin)+roi(1))*pixelsize(1);
% locdat.xnm1=(locs.xpix1(indin)+roi(1))*pixelsize(1);
% % locdat.xnm1=locdat.xnm;
% % w1=1./locdat.xpixerr; w2=1./locdat.xpix2err;
% % locdat.xnm=(locdat.xnm1.*locdat.phot+locdat.xnm2.*locdat.phot2)./(locdat.phot+locdat.phot2);
% % locdat.xnm=(locdat.xnm1.*w1+locdat.xnm2.*w2)./(w1+w2);
% end
% if isfield(locs,'ypix2')
% locdat.ynm2=(locs.ypix2(indin)+roi(2))*pixelsize(end);
% locdat.ynm1=(locs.ypix1(indin)+roi(2))*pixelsize(end);
% % locdat.ynm1=locdat.ynm;
% % w1=1./locdat.ypixerr; w2=1./locdat.ypix2err;
% % locdat.ynm=(locdat.ynm1.*locdat.phot+locdat.ynm2.*locdat.phot2)./(locdat.phot+locdat.phot2);
% % locdat.ynm=(locdat.ynm1.*w1+locdat.ynm2.*w2)./(w1+w2);
% end

fn=fieldnames(locs);
for k=1:length(fn)
    if contains(fn{k},'xpix')
        if contains(fn{k},'err')|| ~strcmp(fn{k}(1),'x')
            offs=0;
        else
            offs=roi(1);
        end  
        newfield=strrep(fn{k},'xpix','xnm');
        locdat.(newfield)=(locs.(fn{k})(indin)+offs)*pixelsize(1);
    end
    if contains(fn{k},'ypix')
        if contains(fn{k},'err') || ~strcmp(fn{k}(1),'y')
            offs=0;
        else
            offs=roi(2);
        end  
        newfield=strrep(fn{k},'ypix','ynm');
        locdat.(newfield)=(locs.(fn{k})(indin)+offs)*pixelsize(1);
    end
end

if isfield(locs,'PSFxpix')
% locdat.PSFxnm=locs.PSFxpix(indin)*pixelsize(1);
else
    locdat.PSFxnm=0*locdat.xnm+100;
end
if isfield(locs,'PSFypix')
%     locdat.PSFynm=locs.PSFypix(indin)*pixelsize(end);
else
    locdat.PSFynm=locdat.PSFxnm;
end
locdat.locprecnm=sqrt((locdat.xnmerr.^2+locdat.ynmerr.^2)/2);
if isfield(locs,'zerr')
    locdat.locprecznm=locs.zerr(indin);
end
if isfield(locs,'znmerr')
    locdat.locprecznm=locs.znmerr(indin);
end

if isfield(locs,'logLikelihood')
    try
    roisize=obj.getPar('loc_ROIsize');
    catch
        roisize=10;
    end
    locdat.LLrel=real(locs.logLikelihood(indin)*2/roisize^2);
end

locdat.filenumber=uint8(0*locdat.xnm+obj.filenumber);
locdat.channel=0*locdat.xnm;
end