function [b,p]=images2beads_globalfit(p)
% determines position of beads and cuts out rois from the bead stacks
% conatainign an individual bead
%initialize some parameters:
fs=p.filtersize; %for smoothing
hfilterim=fspecial('gaussian',2*round(fs*3/2)+1,fs); 
fmax=0;
roisize=p.ROIxy;
roisizeh=min(round(1.5*(p.ROIxy-1)/2),(p.ROIxy-1)/2+1); %create extra space if we need to shift;
rsr=-roisizeh:roisizeh;
filelist=p.filelist;
b=[];
ht=uitab(p.tabgroup,'Title','Files');
tg=uitabgroup(ht);
is4pi=contains(p.modality,'4Pi');

if p.isglobalfit 
   if ischar(p.Tfile)
        l=load(p.Tfile);
        transform=l.transformation;
   else
       transform=p.Tfile;
   end
    p.transformation=transform;
    p.mirror=transform.info{2}.mirror;
else
    p.mirror=false;
end

% load beads in all files
p.multifile=false;
for k=1:length(filelist)
    ax=axes(uitab(tg,'Title',num2str(k)));
%     axis(ax,'image')
    p.fileax(k)=ax;
    if contains(filelist{k},';')
        p.multifile=true;
        ind=strfind(filelist{k},';');
        filelisth=filelist{k}(1:ind-1);
        filelisth2=filelist{k}(ind+1:end);
        [imstack2, p.roi2{k}, p.pixelsize2{k},settings3D2,isem1]=readbeadimages(filelisth2,p);
        [imstack, p.roi{k}, p.pixelsize{k},settings3D,isem2]=readbeadimages(filelisth,p);
        p.emgain=[isem1,isem2];
    else
        filelisth=filelist{k};
        [imstack, p.roi{k}, p.pixelsize{k},settings3D,p.emgain]=readbeadimages(filelisth,p);
        p.roi2=p.roi;
        imstack2=imstack; %for simplicity and to avoid code duplication: when only one channel is considered
    end

    if is4pi
        if ~isempty(settings3D)
            p.settings_3D=settings3D;
            p.settings_3D.file='multifile';
        end
        if isfield(p,'settings_3D') && ~isempty(p.settings_3D) %calibration file: cut out and mirror already here!
            imstack=cutoutchannels(imstack,p.settings_3D);
            imstack2=cutoutchannels(imstack2,p.settings_3D);
            p.roi{k}=[0 0 size(imstack,1) size(imstack,2)]; %check x,y
        else
            disp('no  settings_3D found. Use default. Specify settings file?')
            wx=size(imstack,2)/4;wy=size(imstack,1);
            p.settings_3D=struct('y4pi',[0 0 0 0],'x4pi',[0 wx 2*wx 3*wx], 'width4pi',wx,'height4pi',wy,'mirror4pi',[0 0 0 0],'pixelsize_nm',100,'offset',100,'conversion',0.5);
        end
    end
       
    if isfield(p,'framerangeuse')
        imstack=imstack(:,:,p.framerangeuse(1):p.framerangeuse(end));
        imstack2=imstack2(:,:,p.framerangeuse(1):p.framerangeuse(end));
    end
    
    imstack=imstack-min(imstack(:)); %fast fix for offset;
    imstack2=imstack2-min(imstack2(:)); %fast fix for offset;
    
    %segment beads
    mim=max(imstack,[],3); %maximum intensity projection
    mim=filter2(hfilterim,mim);
    imagesc(ax,mim);
    axis(ax,'image');
    axis(ax,'off')
    title(ax,'Maximum intensity projection')
    if isfield(p,'beadpos') %passed on externally
        maxima=round(p.beadpos{k});
    else
        if isfield(p,'roimask')&&~isempty(p.roimask)
            mim=double(mim).*double(p.roimask);
        end
        maxima=maximumfindcall(mim);
        indgh=maxima(:,1)>p.xrange(1) & maxima(:,1)<p.xrange(2) & maxima(:,2)>p.yrange(1) & maxima(:,2)<p.yrange(2); 
        int=maxima(:,3);
        try
            r1=max(roisize,p.yrange(1)):min(size(mim,1)-roisize,p.yrange(2));
            r2=max(roisize,p.xrange(1)):min(size(mim,2)-roisize,p.xrange(2));
            mimc=mim(r1,r2);
            mmed=myquantile(mimc(:),0.3);
            imt=mimc(mimc<=mmed);
                sm=sort(int);
            mv=mean(sm(end-5:end));
            cutoff=mean(imt(:))+max(2.5*std(imt(:)),(mv-mean(imt(:)))/15);
        catch err
            cutoff=myquantile(mimc(:),.95);
        end
        cutoff=cutoff*p.cutoffrel;
        if any(int>cutoff)
            maxima=maxima(int>cutoff & indgh,:);
        else
            [~,indm]=max(int);
            maxima=maxima(indm,:);
        end
    end
    maximafound=maxima;
    indgoodb=true(size(maxima,1),1);
    
    %remove beads that are closer together than mindistance
    if isfield(p,'mindistance')&&~isempty(p.mindistance)
        for bk=1:size(maxima,1)
            for bl=bk+1:size(maxima,1)
                if  sum((maxima(bk,1:2)-maxima(bl,1:2)).^2)<p.mindistance^2
                    indgoodb(bk)=false;
                    indgoodb(bl)=false;
                end
            end
        end 
       if isfield(p,'settings_3D')
           w=p.settings_3D.width4pi;
           xm=mod(maxima(:,1),w);
           xm(xm>w/2)=xm(xm>w/2)-w;
           indgoodb(abs(xm)<p.mindistance/2)=false;
       end
       hs=size(imstack,1);
       ym=maxima(:,2);
       ym(ym>hs/2)=ym(ym>hs/2)-hs;
       indgoodb(abs(ym)<p.mindistance/2)=false;
       maxima=maxima(indgoodb,:);
    end 
    
  
    if p.isglobalfit
        %calculate in nm on chip (reference for transformation)
        maximanm=(maxima(:,1:2)+p.roi{k}([1 2]));
        if strcmp(transform.unit,'pixel')
            pixfac=[1 1];
            pixfac2=[1 1];
        else
            pixfac=[p.pixelsize{k}(1)*1000 p.pixelsize{k}(end)*1000];
            pixfac2=[p.pixelsize2{k}(1)*1000 p.pixelsize2{k}(end)*1000];
        end
        maximanm(:,1)=maximanm(:,1)*pixfac(1);
        maximanm(:,2)=maximanm(:,2)*pixfac(end);

        %transform reference to target
        indref=transform.getPart(1,maximanm(:,1:2));
        maximaref=maxima(indref,:);
        xy=transform.transformToTarget(2,maximanm(indref,1:2));
        
        %calculate coordinates in pixels on target image
        maximatargetf=[];
        maximatargetf(:,1)=xy(:,1)/pixfac2(1)-p.roi2{k}(1); %now on target chip: use roi2
        maximatargetf(:,2)=xy(:,2)/pixfac2(end)-p.roi2{k}(2);

        maximatar=round(maximatargetf); %round for integer pixel positions
        dxy=maximatargetf-maximatar;   %distance between true transformed bead positions and pixel positions    
    else
        indgoodr=maxima(:,1)>p.xrange(1)&maxima(:,1)<p.xrange(end)&maxima(:,2)>p.yrange(1)&maxima(:,2)<p.yrange(end);
        maxima=maxima(indgoodr,:);
        maximaref=maxima;
        maximatar=maxima;
        dxy=zeros(size(maximatar));
    end
   
    numframes=size(imstack,3);
    bind=length(b)+size(maximaref,1);
    
    %assemble structure with beads
    for l=1:size(maximaref,1)
        b(bind).loc.frames=(1:numframes)';
        b(bind).loc.filenumber=zeros(numframes,1)+k;
        b(bind).filenumber=k;
        b(bind).pos=maximaref(l,1:2);
        b(bind).postar=maximatar(l,1:2);
        b(bind).shiftxy=dxy(l,:);
        try %if transformed beads are outside image, then remove
            b(bind).stack.image=imstack(b(bind).pos(2)+rsr,b(bind).pos(1)+rsr,:);
            b(bind).stack.imagetar=imstack2(b(bind).postar(2)+rsr,b(bind).postar(1)+rsr,:);
            b(bind).stack.framerange=1:numframes;
            b(bind).isstack=true;
        catch err
            b(bind).isstack=false;
        end
        b(bind).roi=p.roi{k};
        b(bind).roi2=p.roi2{k};
        bind=bind-1;
    end
    fmax=max(fmax,numframes);
    hold (ax,'on')
    if p.isglobalfit
        if p.multifile %two channels in different files (e.g. different cameras)
            plot(ax,maxima(:,1),maxima(:,2),'ko',maximafound(:,1),maximafound(:,2),'r.')
            ax2=axes(uitab(tg,'Title',[num2str(k) 't']));
            mim2=max(imstack2,[],3);
            imagesc(ax2,mim2)
            ax2.NextPlot='add';
            axis(ax2,'image');
            axis(ax2,'off')
            plot(ax2,maximatar(:,1),maximatar(:,2),'md',maximafound(:,1),maximafound(:,2),'r.') 
            hold(ax2,'off');
        else %on same camera chip
            plot(ax,maximaref(:,1),maximaref(:,2),'ko',maximatar(:,1),maximatar(:,2),'md',maximafound(:,1),maximafound(:,2),'r.') 
        end
    else
        plot(ax,maxima(:,1),maxima(:,2),'ko',maximafound(:,1),maximafound(:,2),'r.')
    end
    hold (ax,'off')
    drawnow
end
indgoodbead=[b(:).isstack];
b=b(indgoodbead);
p.fminmax=[1 fmax];
p.pathhere=fileparts(filelist{1});
end


