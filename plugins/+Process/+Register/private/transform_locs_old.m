function transform=transform_locs(locData,p)

%get fields

loc.x=locData.getloc('xnm');
loc.y=locData.getloc('ynm');
loc.frame=locData.getloc('frame');
loc.locprec=locData.getloc('locprecnm');
loc.filenumber=locData.getloc('filenumber');
if p.isz
    loc.z=locData.getloc('znm');
else
    loc.psf=locData.getloc('PSFxnm');
end

            
%fitler fields
ofilter=rec.LocalizationFilter;
if p.useroi
    [a,b,ig]=locData.getlocRoi('frame',1);
else
    ig1=ofilter.minMaxFilter(loc.locprec,0,p.maxlocprec);
    if isfield(loc,'z')
        ig2=ofilter.minMaxFilter(loc.z,p.zmin,p.zmax);
    else
        ig2=ofilter.minMaxFilter(loc.psf,0,p.PSFmax);
    end
    ig=ig1&ig2;
end


targetfile=locData.files.file(p.targetselect.value);
reffile=locData.files.file(p.refselect.value);

itarget=loc.filenumber==targetfile.number;
iref=loc.filenumber==reffile.number;

% pixrec=p.pixrec;
pixelsize=reffile.info.pixsize;
%  sizerefpix=[reffile.info.Width reffile.info.Height]*pixelsize*1000;
 roi=reffile.info.roi*pixelsize*1000;
 
 %determine midpoint
 if p.midpcheck
     midp=p.midp*pixelsize*1000;
 else
     switch p.refpart.selection
         case 'all'
             midp=roi(1)+roi(3)/2;
             
         case 'top-bottom'
             midp=roi(2)+roi(4)/2;
         case 'right-left'
             midp=roi(1)+roi(3)/2;
     end
 end
 %
midpoints(1)=roi(1)+roi(3)/2;
midpoints(2)=roi(2)+roi(4)/2;


switch p.targetmirror.selection
    case 'left-right'
        mirror.direction='x';
        
    case 'up-down'
        mirror.direction='y';
    case 'both'
        mirror.direction='xy';
    case 'no mirror'
        mirror.direction='0';
    otherwise
        disp('error')
end
mirror.midpoint(1)=midpoints(1)+roi(3)/4;
mirror.midpoint(2)=midpoints(2)+roi(4)/4;

transform= rec.LocTransform;
transform.midpoint=midp;
transform.refPart=p.refpart.selection;
transform.transformType=p.transform.selection;
transform.transformParameter=p.transformparam;
transform.mirror=mirror;

 
[imr,rangereference]=transform.getSRimage(loc.x,loc.y,ig&iref,'reference',p.pixrec,roi);
[imt,rangetarget]=transform.getSRimage(loc.x,loc.y,ig&itarget,'target',p.pixrec,roi);

% figure(88);
% imagesc(imr);
% figure(89);
% imagesc(imt);

%  p.refpart now different. Calculate targetpart, refpart.
%or: in Transformlocs: set part as property. Then: getPart('reference' etc.)
% in TransfromLocs:
% [iparttarget,mrtarget]=getpart_register(p.targetpart.value,loc.x,loc.y,sizerefpix,midp);
% [ipartref,mrref]=getpart_register(p.refpart.value,loc.x,loc.y,sizerefpix,midp);

% it=ig&itarget&iparttarget;
% ir=ig&iref&ipartref;

% loctarget.x=loc.x(it);
% loctarget.y=loc.y(it);
% loctarget.frame=loc.frame(it);
% loctarget.lp=loc.locprec(it);
% loctref.x=loc.x(ir);
% loctref.y=loc.y(ir);
% loctref.frame=loc.frame(ir);
% loctref.lp=loc.locprec(ir);

%reconstruct SR images

%maybe intransform: getSRimage(loc,indin,'target')
%in there: 
% imr=getsrimage(loctref,p.pixrec,mrref);
% imt=getsrimage(loctarget,p.pixrec,mrtarget);

maxshift=round(p.maxshift/p.pixrec);
ax1=recgui.initaxis(p.resultstabgroup,'shiftcorr');

% s=size(imr)
% ima=zeros(s(1),s(2),3);
% ima(:,:,1)=imr;ima(:,:,2)=imt;
% figure(88);imagesc(ima)
% asdf
[dxpt,dypt]=getShiftCorr(imr,imt,1,maxshift);
dx=dxpt*p.pixrec-rangetarget.x(1)+rangereference.x(1);
dy=dypt*p.pixrec-rangetarget.y(1)+rangereference.y(1);

[loctarget,indt]=transform.getPartLocs(loc,loc.x(ig),loc.y(ig),'target');
[locref,indr]=transform.getPartLocs(loc,loc.x(ig),loc.y(ig),'reference');


% locrefred=reduceFieldLength(locref,ig);
% loctargetred=reduceFieldLength(loctarget,ig);

[iAa,iBa,na]=matchlocsall(locref,loctarget,dx,dy,p.maxshiftregister,p.maxlocsused);

% transform.findTransform(locref.x(iAa),locref.y(iAa),loctarget.x(iBa),loctarget.y(iBa))
transform.findTransform(locref.x(iAa),locref.y(iAa),loctarget.x(iBa),loctarget.y(iBa))

if p.showresults
    recgui.initaxis(p.resultstabgroup,'scatter')
    [xa, ya]=transform.transformCoordinates((loctarget.x(iBa)),(loctarget.y(iBa)),'target');
%  [xa, ya]=transform.transformCoordinates((loc.x(indr(iBa))),(loc.y(indt(iBa))),'target');
 
   dx=xa-locref.x(iAa);
   dy=ya-locref.y(iAa);
   
%    dx=loctarget.x(iBa)-locref.x(iAa);
%    dx=loctarget.y(iBa)-locref.y(iAa);
%    figure(88)
   dscatter(dx,dy)
   title({['number of anchor points: ' num2str(length(iBa)) ' of ' num2str(length(iAa)+length(na))],['dx= ' num2str(std(dx),3) ' nm, dy= ' num2str(std(dy),3) ' nm']});
 ax3=recgui.initaxis(p.resultstabgroup,'hist');
 hist(dx,50)
end



% [f,path]=uiputfile(p.Tfile);
% if f
%output: transform object
%     transform.T=mytform;
%     transform.midpoint=midp;
%     transform.refpart=p.refpart.value;
%     transform.refmirror=0; %%not implemented
%     transform.targetpart=p.targetpart.value;
%     transform.targetmirror=0; %%not implemented
%     save([path filesep f],'transform')
% end
% 
% x=loc.x(ig);
% y=loc.y(ig);
% f=loc.frame(ig);
% fn=loc.filenumber(ig);


% mx=myquantile(x,0.9




function outim=getsrimage(l,pixrec,m)

[outim,dxo,dyo]=myhist2(l.x,l.y,pixrec,pixrec,m.x,m.y);
   

