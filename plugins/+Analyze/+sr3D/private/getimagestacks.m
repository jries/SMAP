function [beads,allstacks]=getimagestacks(obj,p,beads)

roisize=p.roisize;
% halfroisize=(roisize-1)/2;
roifac=1.5;
roisizebig=(roifac*roisize);
halfroisizebig=round((roisizebig-1)/2); %room for shifting
roisizebig=2*halfroisizebig+1;

f0=[beads(:).f0];
f0r=round(f0);

fminmax=[min(obj.locData.loc.frame) max(obj.locData.loc.frame)];
if p.beaddistribution.Value==1&&~p.alignz
    storeframes=p.roiframes;
else
    storeframes=p.roiframes+300/p.dz;
end
halfstoreframes=round((storeframes-1)/2);
storeframes=halfstoreframes*2+1;

if p.beaddistribution.Value==2|| (p.beaddistribution.Value==1&&~p.alignz)
    fzero=halfstoreframes+1;
else%if on glass
    fzero=round(median(f0(~isnan(f0))));
end

framerange0=max(fminmax(1),fzero-halfstoreframes):min(fzero+halfstoreframes+1,fminmax(2));

files=obj.locData.files.file;
indbad=false(length(beads),1);
allstacks=[];
 for k=1:length(files) 
    if isfield(files(k).info,'imagefile')
        filename=files(k).info.imagefile;
    else
        filename=files(k).info.basefile;
    end
    roi=files(k).info.roi;
    campix=files(k).info.cam_pixelsize_um;
    il=getimageloader(obj,filename);
    imstackadu=il.getmanyimages([],'mat');
    il.close;
                
    imstack=(double(imstackadu)-il.metadata.offset)*il.metadata.pix2phot;
    si=size(imstack);
    if isempty(allstacks)
        allstacks=zeros(si(1),si(2),si(3),length(files));
    end
    allstacks(1:si(1),1:si(2),1:si(3),k)=imstack;
    sim=size(imstack);
    thisfile=find([beads(:).filenumber]==k);
%                 f=figure(88);
%             hold off    
%             imagesc(imstack(:,:,40));
            
    for s=length(thisfile):-1:1
        beadnumber=thisfile(s);
        pos=round(beads(beadnumber).pos(1:2)./campix/1000-roi(1:2));
        if pos(1)>halfroisizebig&&pos(1)<sim(2)-halfroisizebig&&pos(2)>halfroisizebig&&pos(2)<sim(1)-halfroisizebig
            if p.beaddistribution.Value==1&&~p.alignz %on glass
                frange=framerange0;
            else
                frange=max(fminmax(1),f0r(beadnumber)-halfstoreframes):min(fminmax(2),f0r(beadnumber)+halfstoreframes);
            end
            if length(frange)<halfroisizebig
                indbad(beadnumber)=true;
                continue
            end
            smallfr=imstack(pos(2)-halfroisizebig:pos(2)+halfroisizebig,pos(1)-halfroisizebig:pos(1)+halfroisizebig,frange);
            smallf=rescalestack(smallfr);
            stack=zeros(roisizebig,roisizebig,storeframes)+NaN;
            % XXXXX gel: or always: center fzero or similar.
            if frange(1)==fminmax(1) && frange(end)~=fminmax(2)
                stack(:,:,end-size(smallf,3)+1:end)=smallf;
            else
                stack(:,:,1:size(smallf,3))=smallf;
            end
            beads(beadnumber).stack.image=stack;
            beads(beadnumber).stack.framerange=frange;
            
%             figure(88)
% hold on
% plot(pos(1),pos(2),'wo')
%             figure(89)
%             imagesc(smallf(:,:,12));
%             title(pos)
            
            
%             imageslicer(stack,'Parent',f)
            
        else
            indbad(beadnumber)=true;
        end
                    
    end
 end
beads(indbad)=[];

end

function out=rescalestack(in)
% out=in;
% return;

s=size(in);
midp=round((s(1)-1)/2+1);
range=3;
center=squeeze(sum(sum(in(midp-range:midp+range,midp-range:midp+range,:))));

% out=in/max(in(:));
try
x=(1:s(3))';
fp=fit(x,center,'gauss1','Lower',[0 1 8],'Upper',[inf length(x) length(x)*2]);
int=fp(x);
corr=int./center;
for k=s(3):-1:1
    out(:,:,k)=in(:,:,k)*corr(k);
end
catch err
    err
    out=in;
end
% center2=squeeze(sum(sum(out(midp-1:midp+1,midp-1:midp+1,:))));
end
