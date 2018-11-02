function [locsout,possites]=simulatelocs2(p, colour)
           [poslabels,possites]=getlabels(p, colour);
%            p.maxframe=100000;
           posreappear=reappear(poslabels,p.blinks(colour),p.maxframes(colour));
           
%            p.lifetime=2;
           photonsperframe=p.photons(colour)/p.lifetime(colour);
           posphot=getphotons(posreappear,photonsperframe,p.lifetime(colour));
           
           locs=locsfrompos(posphot,p, colour);
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
timeonall=[lenfirst; timeon];
photons=timeonall*photonsperframe;
photonsr=poissrnd(photons);
locso.phot=photonsr;
locso.frame=[locs.frame; frame];

end

function [locs,possites]=getlabels(p, colour)
%gets position of labels either from image or from coordinates
% fields p. :
% coordinatefile, se_sitefov, numberofsites(x,y), labeling_efficiency, randomxy,
% randomxyd

if p.usePes
    switch colour
        case 1  
            [image, nm] = p.obj.getPar('proteinES').(p.selectPro.selection).getTimeImg(p.time, p.viewType.selection, p.tif_imagesize);
            image = imrotate(image, 180);
            p.Mean = nm;
            p.Std = 0;
        case 2
            [image, nm] = p.obj.getPar('proteinES').(p.pro2nd.selection).getTimeImg(p.time, p.viewType.selection, p.tif_imagesize);
            image = imrotate(image, 180);
            p.Mean = nm;
            p.Std = 0;
        otherwise
    end
else
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
            p.Mean = 150;
            p.Std = 20;

        case '.mat'
            l=load(p.coordinatefile);

            if isfield(l,'image')
                image=l.image;
            else
                locsall=copyfields([],l,{'x','y','z'});
            end
        case '.m'
            cf=pwd;
            [ph,fh]=fileparts(p.coordinatefile);
            cd(ph)
            l=eval(fh);
            cd(cf);
    %         l=eval(p.coordinatefile);
            if isfield(l,'image')
                image=l.image;
            else
                locsall=copyfields([],l,{'x','y','z','channel'}); % channel is added
            end
    % add a new way to get localizations
        case {'.loc'}
            image=imread(p.coordinatefile);
            img=sum(image,3)/size(image,3); %binarize
            image=double(img)/255;
        otherwise
            display('file not identified selected')
            return
    end % after this block, you will get either image or locsall, depending on the type of your input
    %if ~isfield(locsall,'z')
    %    locsall.z=0*locsall.x;
    %end    
end


if isempty(p.se_sitefov)
    p.se_sitefov=500;
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
for k=numberofsites:-1:1
    xh=mod(k-1,numberofrows)+1;
    yh=ceil(k/numberofrows);
    if ~isempty(image)
        locsh=locsfromimage(image,p,colour);
    else
        locsh=labelremove(locsall,p.labeling_efficiency); % give all the coordinates of I dots and the lableling efficiency (P(ref)), and then randomly generating I even probabilities (P(gi)), keep P(gi) if the P(gi) <= P(ref)
    end

    numlocs=length(locsh.x);
    locsh.x=reshape(locsh.x,numlocs,1);
    locsh.y=reshape(locsh.y,numlocs,1);
    locsh.z=reshape(locsh.z,numlocs,1);
    locsh.channel=reshape(locsh.channel,numlocs,1); %added
    if p.randomrot
        angle=2*pi*rand(1);
        [locsh.x,locsh.y]=rotcoord(locsh.x,locsh.y,angle);
    else
        angle=0;
    end
    
    if p.randomxy 
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
%     figure(89);plot(locs(k).x,locs(k).y,'*')
%     waitforbuttonpress
    
end
end


function locs=locsfromimage(image,p, colour)
    % adapted from applyFilter
    
    finalNumMol = round(random('norm', p.Mean, p.Std));
   
    norFilter = image/max(image(:)); % = normalized filter. the maximun of a pixel value is set 1
    
    ExpP = numel(image); % number of image pixels
    ActP = sum(norFilter(:)); % sum of all pixel values in the normalized filter
    numInput = (finalNumMol*ExpP)/ActP;
    numInput = ceil(numInput * 2);
    
    setRange = length(image);
    x = rand(1,numInput)*setRange;
    y = rand(1,numInput)*setRange;
    pr = rand(1,numInput);
    
    idx = sub2ind(size(norFilter), ceil(x), ceil(y));
    
    filteredIdx = pr <= norFilter(idx);
    
    xy = [x; y];
    Filtered = xy(:,filteredIdx);
    Filtered = Filtered(:,1:finalNumMol);
    locs = [];
    locs.x = Filtered(1,:)'-setRange/2;
    locs.y = Filtered(2,:)'-setRange/2;
    locs.z = repelem(0,length(Filtered));
    locs.channel = repelem(colour,length(Filtered));
end

function locs=labelremove(locin,p)
numl=length(locin.x);
indin=rand(numl,1)<=p;
locs=copystructReduce(locin,indin);
end

function locso=reappear(locs,numblinks,maxframe)
for k=length(locs):-1:1
    locso(k)=reappeari(locs(k),numblinks,maxframe);
end
end


function locso=reappeari(locs,numblinks,maxframe)
fn=fieldnames(locs);
numlocs=length(locs.(fn{1}));
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

% if nargin>2
locso.frame=double(ceil(rand(length(indout),1)*maxframe));
% end
end


function locs=locsfrompos(locsi,p, colour)
for k=length(locsi):-1:1
    locs(k)=locsfromposi(locsi(k),p, colour);
end
end

function locs=locsfromposi(locsi,p, colour)
    numlocs=length(locsi.x);
%     phot=exprnd(p.photons,numlocs,1);
    phot=locsi.phot;
    
    a=100;
    PSF=100;
    sa=PSF+a/12;
    phot(phot<10)=10;
    indin=phot>=10;
    numlocs=sum(indin);
    %MOrtensen
    locprecnm=sqrt(sa^2./phot.*(16/9+8*pi*sa^2*p.background(colour)^2./phot/a^2));
    locs.phot=single(phot(indin));
    locs.bg=single(locprecnm(indin)*0+p.background(colour));
    locs.locprecnm=single(locprecnm(indin));
%     locs.frame=double(ceil(rand(numlocs,1)*p.maxframe));
    locs.xnm=single(locsi.x(indin)+randn(numlocs,1).*locprecnm(indin));
    locs.ynm=single(locsi.y(indin)+randn(numlocs,1).*locprecnm(indin));
    locs.znm=single(locsi.z(indin)+randn(numlocs,1).*locprecnm(indin)*3);
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