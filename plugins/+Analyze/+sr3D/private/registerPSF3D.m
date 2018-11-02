function [imout,shiftedstackn,shift,indgood]=registerPSF3D(imin,p,axs)
if nargin<3
    axs={};
end
% perform correlation on:
% p.xrange
% p.yrange
% p.framerange
% can do 2D if length(p.framerange)=1
%cutout small volumes
if ~isfield(p,'xrange')
    p.xrange=1:size(imin,1);
end
if ~isfield(p,'yrange')
    p.yrange=1:size(imin,2);
end
if ~isfield(p,'framerange')
    p.framerange=1:size(imin,3);
end

numbeads=size(imin,4);
if numbeads==1
    imout=imin;
    shiftedstackn=imin;
    shift=[0 0 0 ];
    indgood=true;
    return
end
if ~p.beadfilterf0
    p.framerange=3:size(imin,3)-3;
end
smallim=zeros(length(p.xrange),length(p.yrange),length(p.framerange),size(imin,4));

if ~isempty(p.zshiftf0)
    zshiftf0=p.zshiftf0;
else
    zshiftf0=zeros(numbeads,1);
end

for k=1:size(imin,4)
    frh=round(p.framerange-zshiftf0(k));
    smallim(:,:,:,k)=imin(p.xrange,p.yrange,frh,k);
end
avim=nanmean(imin,4);

xn=1:size(imin,1);yn=1:size(imin,2);zn=1:size(imin,3);
[Xq,Yq,Zq]=meshgrid(yn,xn,zn);
% meanim=zeros(size(Xq));
meanim=[];
refim=avim(p.xrange,p.yrange,p.framerange);

simin=size(imin);
shiftedstack=zeros(simin(1),simin(2),simin(3),numbeads)+NaN;

for k=1:numbeads
%     shift(k,:)=get3Dcorrshift(avim,smallim(:,:,:,k));
    goodframes=squeeze(nansum(nansum(smallim(:,:,:,k),1),2))>0;
    if p.alignz
        [shift(k,:),cc(k)]=get3Dcorrshift(refim(:,:,goodframes),smallim(:,:,goodframes,k));
    else
        if any(goodframes)
            [shift(k,:),cc(k)]=get2Dcorrshift(refim(:,:,goodframes),smallim(:,:,goodframes,k));
        else
            shift(k,:)=[0 0 0];cc(k)=NaN;
        end
    end
    
    shiftedh=interp3(imin(:,:,:,k),Xq-shift(k,2),Yq-shift(k,1),Zq-shift(k,3)-double(zshiftf0(k)),'cubic',0);
%     if ~isempty(meanim)
%         ratio=shiftedh(p.xrange,p.yrange,p.framerange)./meanim(p.xrange,p.yrange,p.framerange);
%         factor=nanmedian(ratio(:));
%     else
%         factor=nanmax(shiftedh(:));
%     end
%     if factor==0
%         factor=NaN;
%     end
    shiftedstack(:,:,:,k)=shiftedh;%/factor;
%     shiftedh(isnan(shiftedh))=0;
%     meanim=shiftedh+meanim;
    meanim=nanmean(shiftedstack(:,:,:,1:k),4);
    meanim(isnan(meanim))=avim(isnan(meanim));
    
    refim=meanim(p.xrange,p.yrange,p.framerange);
end

shiftedstackn=normalizstack(shiftedstack,p);


[indgood]=getoverlap(shiftedstackn,shift,p,true(1,size(shiftedstackn,4)));
[indgood]=getoverlap(shiftedstackn,shift,p,indgood);
[indgood,res,normglobal,co,cc2]=getoverlap(shiftedstackn,shift,p,indgood);

shiftedstackn=shiftedstackn/normglobal;

imout=nanmean(shiftedstackn(:,:,:,indgood),4);
shiftedstackn(1,end,:,~indgood)=nanmax(shiftedstackn(:));
shiftedstackn(1,:,1,~indgood)=nanmax(shiftedstackn(:));

if length(axs)>0
hold(axs{1},'off')
plot(axs{1},(res(~indgood)),cc2(~indgood),'g*');title(co)
hold(axs{1},'on')
plot(axs{1},(res(indgood)),cc2(indgood),'rx');title(co)
xlabel(axs{1},'residulas')
ylabel(axs{1},'cross-correlation value')
end

if length(axs)>1
    imageslicer(vertcat(avim,imout),'Parent',axs{2}.Parent)
end
%later: include 2x upscaling here

%  shift=-shift;
% for k=numbeads:-1:1
%     
%     
% end
% imout=meanim/numbeads;

% f=figure(99);
% delete(f.Children)
% % avim(2,2,2)=max(avim(:));
% imageslicer(vertcat(avim,imout,shiftedstack(:,:,:,1),imin(:,:,:,1)),'Parent',f);
% imageslicer(shiftedstack);
% average
%shift of everything to average
%shift volumea


end

function [indgood,res,normamp,co,cc]=getoverlap(shiftedstackn,shift,p,indgood)
%re-calculate CC! bisas for first.
refimn=nanmean(shiftedstackn(p.xrange,p.yrange,p.framerange,indgood),4);
% refimb=nanmean(shiftedstackn(:,:,:,indgood),4);

for k=size(shiftedstackn,4):-1:1
    imh=shiftedstackn(p.xrange,p.yrange,p.framerange,k);
    badind=isnan(imh)|isnan(refimn);
    cc(k)=sum(refimn(~badind).*imh(~badind))/(sum(refimn(~badind))*sum(imh(~badind)))*sum(~badind(:));
%     imhb=shiftedstackn(:,:,:,k);
%     badindb=isnan(imhb)|isnan(refimb);
%     ccb(k)=sum(refimb(~badindb).*imhb(~badindb))/(sum(refimb(~badindb))*sum(imhb(~badindb)))*sum(~badindb(:));
%     
end

normamp=nanmax(refimn(:));
shiftedstackn=shiftedstackn/normamp;
refimn=refimn/normamp;
for k=size(shiftedstackn,4):-1:1
     sim=shiftedstackn(p.xrange(2:end-1),p.yrange(2:end-1),p.framerange,k);
     dv=(refimn(2:end-1,2:end-1,:)-sim).^2;
    res(k)=sqrt(nanmean(dv(:)));

    
%          simb=shiftedstackn(:,:,:,k);
%      dvb=(refimb(:,:,:)-simb).^2;
%     resb(k)=sqrt(nanmean(dvb(:)));
end
rescc=res./cc;
rescc(abs(shift(:,1))>3|abs(shift(:,2))>3)=NaN;
[a,b]=robustMean(rescc(cc>0));
if isnan(b)
    a=nanmean(rescc);b=nanstd(rescc);
end
co=a+3.*b;
indgood=rescc<co;
end

function out=normalizstack(in,p)
sin=size(in);
  midp=round((length(p.xrange)+1)/2);
    xr=p.xrange(midp-3:midp+3);yr=p.yrange(midp-3:midp+3);
if p.beadfilterf0 
    out=0*in+NaN;
  
    for k=1:sin(4)
        imh=in(xr,yr,p.framerange,k);
        nm=nanmean(imh(:));
        if nm>0
        out(:,:,:,k)=in(:,:,:,k)/nm;
        end
    end
else %use fitting
%     x0=ones(1,size(in,4));
%     fitp=lsqnonlin(@alignstacks,x0,0*x0,[],[],in)

    inh=in;
    out=0*in+NaN;
%     xr=p.xrange(2:end-1);yr=p.yrange(2:end-1);
    for iter=1:4
        meanim=nanmean(inh,4);
        for k=1:sin(4)
            imh=inh(:,:,:,k);
            if all(isnan(imh))
                continue
            end
            ims=imh(xr,yr,p.framerange);
            meanims=meanim(xr,yr,p.framerange);
            isn=isnan(ims)|isnan(meanims);
            intcutoff=meanims>myquantile(meanims(:),0.75);
            indg=~isn&intcutoff;
            ratio=ims(indg)./meanims(indg);
            
            factor=nanmedian(ratio(:));
            factor2=sum(ims(indg))/sum(meanims(indg));
            disp([factor  factor2])
            if factor>0
                out(:,:,:,k)=imh/factor;
            end
        end
        inh=out;
    end
%         iii=squeeze(out(round(sin(1)/2),round(sin(1)/2),:,:));
%         iii2=squeeze(in(round(sin(1)/2),round(sin(1)/2),:,:));
%     figure(88);subplot(2,1,1);plot(iii);subplot(2,1,2);plot(iii2);hold off
end

end

