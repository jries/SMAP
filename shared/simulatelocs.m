function [locsout,possites,parameters]=simulatelocs(p, colour)
        if nargin<2
            colour=1;
        end
           [poslabels,possites,parameters]=getlabels(p, colour);
           
           if ~isfield(p,'model')
               p.model.selection='simple';
           end
           if ~isfield(p,'blinks')
               p.blinks=0;
           end
           if ~isfield(p,'maxframes')
               p.maxframes=100;
           end
           if ~isfield(p,'lifetime')
               p.lifetime=1;
           end
           if ~isfield(p,'background')
               p.background=0;
           end
           posreappear=getblinks(poslabels,p.model.selection,p.blinks,p.maxframes);
           
%            p.lifetime=2;
           photonsperframe=p.photons/p.lifetime;
           posphot=getphotons(posreappear,photonsperframe,p.lifetime);
           
           locs=locsfrompos(posphot,p);
           locsout=singlelocs(locs);
end

function posphot=getphotons(locs,photonsperframe,lifetime)
for k=length(locs):-1:1
    posphot(k)=getphotonsi(locs(k),photonsperframe,lifetime);
end
end

function locso=getphotonsi(locs,photonsperframe,lifetime)
%calculate intensities:
%total lifetime
fn=fieldnames(locs);
numlocs=length(locs.(fn{1}));
lifetime=exprnd(lifetime,numlocs,1);
lenfirst=rand(numlocs,1);
lenframe=ceil(lifetime-lenfirst);

indextra=zeros(sum(lenframe),1);
frame=zeros(sum(lenframe),1);
timeon=zeros(sum(lenframe),1);


idx=1;

for k=1:length(lenframe)
    indextra(idx:idx+lenframe(k)-1)=k;
%     frame(idx:idx+lenframe(k)-1)=locs.frame(k)+(1:lenframe(k));
    T=lifetime(k)-lenfirst(k);
%     while (T>1)
%         timeon(
%         T=T-1;
%     end
    
    for l=1:lenframe(k)
        frame(idx+l-1)=locs.frame(k)+l;
        if T>1
            timeon(idx+l-1)=1;
            T=T-1;
        else
            timeon(idx+l-1)=T;
        end
    end
    idx=idx+lenframe(k);
end



indout=[(1:numlocs)'; indextra];
for k=1:length(fn)
    locso.(fn{k})=locs.(fn{k})(indout);
end
timeonall=[min(lenfirst,lifetime); timeon];
photons=timeonall*photonsperframe;
photonsr=poissrnd(photons);
locso.phot=photonsr;
locso.frame=[locs.frame; frame];
%remove double locs in one frame coming from same fluorophore
A=horzcat(locso.fluorophore, locso.frame);
[~,ia,ic]=unique(A,'rows');
ind=false(length(ia),1);
ind(ia)=true;
locso=copystructReduce(locso,ind);

end

function [locs,possites,parameters]=getlabels(p, colour)
%gets position of labels either from image or from coordinates
% fields p. :
% coordinatefile, se_sitefov, numberofsites(x,y), labeling_efficiency, randomxy,
% randomxyd
if ischar(p.coordinatefile)
paths =  strsplit(p.coordinatefile, '|');
switch colour
    case 1
        p.coordinatefile = paths{1};
    case 2
        p.coordinatefile = paths{2};
    otherwise
end
        
[~,~,ext]=fileparts(p.coordinatefile);% Get extension of the specified file
image=[];
locsall=[];
locfun=[];
switch ext
    case {'.txt','.csv'}
        plocs=readtable(p.coordinatefile);
        plocsa=table2array(plocs);
        locsall.x=plocsa(:,1);
        locsall.y=plocsa(:,2);
        if size(plocsa,2)>2
            locsall.z=plocsa(:,3);
        end
        locsall=copyfields(locsall,plocs,{'x','y','z'});
    case {'.tif','.png'}
%         locs=getlabelstiff(obj,p);
        image=imread(p.coordinatefile);
        img=sum(image,3)/size(image,3); %binarize
        image=double(img)/255;
        
    case '.mat'
        l=load(p.coordinatefile);
        
        if isfield(l,'image')
            image=l.image;
        else
            locsall=copyfields([],l,{'x','y','z','channel'});
        end
    case '.m'
            locfun=p.coordinatefile;
% add a new way to get localizations
    case {'.loc'}
        image=imread(p.coordinatefile);
        img=sum(image,3)/size(image,3); %binarize
        image=double(img)/255;
    otherwise
        display('file not identified selected')
        return
end % after this block, you will get either image or locsall, depending on the type of your input
else %pass on matlab variable
    if isstruct(p.coordinatefile)
        locsall=copyfields([],p.coordinatefile,{'x','y','z','channel'});
    else
        image=p.coordinatefile;
    end
        
end


if ~isfield(p,'se_sitefov') || isempty(p.se_sitefov)
    if isfield(p,'size')
        p.se_sitefov=p.size;
    else
        p.se_sitefov=500;
    end
end

if ~isfield(p,'labeling_efficiency') || isempty(p.labeling_efficiency)
    p.labeling_efficiency=1;
end

distsites=p.se_sitefov;

if length(p.numberofsites)>1
    numberofrows=p.numberofsites(2);
    numberofsites=p.numberofsites(1)*p.numberofsites(2);
else
    numberofrows=ceil(32000/p.se_sitefov);
    numberofsites=p.numberofsites;
end

% numeroflines=ceil(p.numberofsites/numberofrows);
fieldstosave={'labeling_efficiency','model','blinks','lifetime','photons','background','maxframes','coordinatefile'};
psave=copyfields([],p,fieldstosave);
for k=numberofsites:-1:1
    xh=mod(k-1,numberofrows)+1;
    yh=ceil(k/numberofrows);
    
    phere=psave;
    if ~isempty(image)
        locsh=locsfromimage(image,p);
        
    elseif ~isempty(locfun)
        [locsd,ph]=getcoordm(p);
        locsh=labelremove(locsd,p.labeling_efficiency);
        phere.model=ph;
    else
        locsh=labelremove(locsall,p.labeling_efficiency); % give all the coordinates of I dots and the lableling efficiency (P(ref)), and then randomly generating I even probabilities (P(gi)), keep P(gi) if the P(gi) <= P(ref)
    end
    if ~isfield(locsh,'z')
        locsh.z=0*locsh.x;
    end
    numlocs=length(locsh.x);
    locsh.x=reshape(locsh.x,numlocs,1);
    locsh.y=reshape(locsh.y,numlocs,1);
    locsh.z=reshape(locsh.z,numlocs,1);
    locsh.channel=reshape(locsh.channel,numlocs,1); %added
    if isfield(p,'randomrot') && p.randomrot
        angle=2*pi*rand(1);
        phere.angle=angle;
        [locsh.x,locsh.y]=rotcoord(locsh.x,locsh.y,angle);
        if isfield(locsh,'z')
            dtheta=p.randomrotangle/360*2*pi;
            thetas=dtheta*rand(1);
            [locsh.x,locsh.z]=rotcoord(locsh.x,locsh.z,thetas);
            angle=2*pi*rand(1);
            [locsh.x,locsh.y]=rotcoord(locsh.x,locsh.y,angle);
            phere.theta=thetas;
        end
    else
        angle=0;
    end
    
    if isfield(p,'randomxy') && p.randomxy 
        if length(p.randomxyd)==1
            zspread=p.randomxyd;
        else
            zspread=p.randomdxyd(2);
            p.randomdxyd=p.randomdxyd(1);
        end
        dx=(rand(1)-.5)*p.randomxyd*2;dy=(rand(1)-.5)*p.randomxyd*2;dz=(rand(1)-.5)*zspread*2;
        locsh.x=locsh.x+dx;
        locsh.y=locsh.y+dy;
        locsh.z=locsh.z+dz;
        phere.dxyz=[dx dy dz];
    else
        dx=0;
        dy=0;
        dz=0;
    end
    
    locs(k).x=locsh.x+xh*distsites;
    locs(k).y=locsh.y+yh*distsites;
    locs(k).z=locsh.z;
    locs(k).channel=locsh.channel;% added
    locs(k).angle=angle*ones(size(locsh.x));
    locs(k).dx_gt=dx*ones(size(locsh.x));
    locs(k).dy_gt=dy*ones(size(locsh.x));
    locs(k).dz_gt=dz*ones(size(locsh.x));
    
    possites(k).x=xh*distsites;
    possites(k).y=yh*distsites;
    parameters(k)=phere;
%     figure(89);plot(locs(k).x,locs(k).y,'*')
%     waitforbuttonpress
    
end
end


function locs=locsfromimage(image,p)

if p.tif_numbermode.Value==1
    density=p.tif_density;
else
    %calcualte density from number of locs
%     pdensity=mean(image(:))/(p.tif_imagesize/1000)^2;
    density=p.tif_density/mean(image(:))/(p.tif_imagesize/1000)^2;
end
numtot=round(density*(p.tif_imagesize/1000)^2);

x=rand(numtot,1);
y=rand(numtot,1);


xpix=ceil(x*size(image,1));
ypix=ceil(y*size(image,2));
linind=sub2ind(size(image),xpix,ypix);
keep=image(linind)>rand(numtot,1);

locs.x=(x(keep)-0.5)*p.tif_imagesize;
locs.y=(y(keep)-0.5)*p.tif_imagesize;
% pixf=size(image,1)/(p.tif_imagesize/1000);
% densitypixel=density/pixf^2;

end

function locs=labelremove(locin,p)
numl=length(locin.x);
indin=rand(numl,1)<=p;
locs=copystructReduce(locin,indin);
end

function locso=getblinks(locs,model,numblinks,maxframe)
lo=0;
ln=0;
for k=length(locs):-1:1
    locso(k)=getblinksi(locs(k),model,numblinks,maxframe);
    lo=lo+length(locs(k).x);
    ln=ln+length(locso(k).x);
end

 disp(['occurance per fluorophore: ' num2str(ln/lo)])
end


function locso=getblinksi(locs,model,numblinks,maxframe)
fn=fieldnames(locs);
numlocs=length(locs.(fn{1}));
switch model
    case {'simple','PAFP'}
        locs.frame=double(ceil(rand(numlocs,1)*maxframe));
        fn=fieldnames(locs);
        numn=round(exprnd(numblinks,numlocs,1));
        indextra=zeros(sum(numn),1);
        idx=1;
        for k=1:length(numn)
            indextra(idx:idx+numn(k)-1)=k;
            idx=idx+numn(k);
        end
        indout=[(1:numlocs)'; indextra];
        for k=1:length(fn)
            locso.(fn{k})=locs.(fn{k})(indout);
        end
        
        if contains(model,'simple')
            locso.frame=double(ceil(rand(length(indout),1)*maxframe));
%             locso.fluorophore=indout;
        else % 'PAFP'
            frames=locs.frame(1:numlocs);
            for k=numlocs+1:length(indout)
                timescale=10;
                df=ceil(exprnd(timescale));
                frames(indout(k))=frames(indout(k))+df;
                locso.frame(k)=frames(indout(k));
            end
%             locso.fluorophore=indout;
            
        end
        
    case 'Dye'
        timescale=maxframe/(numblinks+1)/2;
        frame=(exprnd(timescale,numlocs,1));
        on=true(size(frame));
        onind=find(on);
        indout=(1:numlocs)';
        pbl=1/(numblinks+1);
        framenow=frame;
        numon=sum(on);
        while numon>length(on)/1000 %.1% left
            %bleach
            indbl=rand(numon,1)<pbl;
            on(onind(indbl))=false;
            onind=find(on);
            numon=sum(on);
            indout=[indout;onind];
            framenow(on)=framenow(on)+(exprnd(timescale,numon,1));
            frame=[frame; framenow(on)];   

        end
        for k=1:length(fn)
            locso.(fn{k})=locs.(fn{k})(indout);
        end
%         locso.fluorophore=indout;
        
        %distribute equal, stretch to maxframes
        [~,sortind]=sort(frame);
        locso.frame=round(sortind/length(sortind)*maxframe);   
    otherwise 
        disp('model not implemented')
end
locso.fluorophore=indout;
% end
end


function locs=locsfrompos(locsi,p)
for k=length(locsi):-1:1
    locs(k)=locsfromposi(locsi(k),p);
end
end

function locs=locsfromposi(locsi,p)
%     numlocs=length(locsi.x);
%     phot=exprnd(p.photons,numlocs,1);
    phot=locsi.phot;
    
    a=100;
    PSF=100;
    zfactor=1;
%     sa=PSF+a/12;
    phot(phot<10)=10;
    indin=phot>=10;
    numlocs=sum(indin);
    %MOrtensen
%     p.background: photons, b^2=p.background: variance in Mortenson
    locpthompson=sqrt((PSF^2+a^2/12)./(phot)+8*pi*(PSF^2).^2.* p.background./(phot).^2/a^2);
    locprecnm=sqrt((PSF^2+a^2/12)./phot.*(16/9+8*pi*(PSF^2+a^2/12)*p.background./phot/a^2));
    
    
    locs.phot=single(phot(indin));
    locs.bg=single(locprecnm(indin)*0+p.background);
    locs.locprecnm=single(locprecnm(indin));
%     locs.frame=double(ceil(rand(numlocs,1)*p.maxframe));
    locs.xnm=single(locsi.x(indin)+randn(numlocs,1).*locprecnm(indin));
    locs.ynm=single(locsi.y(indin)+randn(numlocs,1).*locprecnm(indin));
    locs.znm=single(locsi.z(indin)+randn(numlocs,1).*locprecnm(indin)*zfactor);
    locs.xnm_gt=single(locsi.x(indin));
    locs.ynm_gt=single(locsi.y(indin));
    locs.znm_gt=single(locsi.z(indin));
    locs.dxnm_gt=single(locsi.dx_gt(indin));
    locs.dynm_gt=single(locsi.dy_gt(indin));
    locs.dznm_gt=single(locsi.dz_gt(indin));
    locs.angle=single(locsi.angle(indin));
    locs.frame=single(locsi.frame(indin));
    locs.channel=single(locsi.channel(indin));  %added
end


function locso=singlelocs(locs)
numl=0;
fn=fieldnames(locs(1));
for k=1:length(locs)
    numl=numl+length(locs(k).(fn{1}));
end

% locso=locs(1);
% locso.(fn{1})(1)=0;
for k=1:length(fn)
    locso.(fn{k})(numl,1)=locs(1).(fn{k})(1); %initialize
end

ind=1;
for l=1:length(locs)
    numlh=length(locs(l).(fn{1}));
    for k=1:length(fn)
        locso.(fn{k})(ind:ind+numlh-1)=locs(l).(fn{k});
    end
    ind=ind+numlh;
end
end

function [out,parameters]=getcoordm(p)
        cf=pwd;
        [ph,fh]=fileparts(p.coordinatefile);
        cd(ph)
        try
        [l,parameters]=eval(fh);
        catch
            [l]=eval(fh);
            parameters=[];
        end
        cd(cf);
%         l=eval(p.coordinatefile);
        if isfield(l,'image')
           
            out=locsfromimage(l.image,p);
        else
            out=copyfields([],l,{'x','y','z','channel'}); % channel is added
        end
end