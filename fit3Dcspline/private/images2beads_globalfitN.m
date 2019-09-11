function [b,p]=images2beads_globalfitN(p)
% addpath('bfmatlab')
fs=p.filtersize;
h=fspecial('gaussian',2*round(fs*3/2)+1,fs);
fmax=0;
roisize=p.ROIxy;
roisizeh=round(1.5*(p.ROIxy-1)/2); %create extra space if we need to shift;
rsr=-roisizeh:roisizeh;
filelist=p.filelist;
b=[];
ht=uitab(p.tabgroup,'Title','Files');
tg=uitabgroup(ht);

% transform=p.transformation;

% offset=3;
for k=1:length(filelist)
    ax=axes(uitab(tg,'Title',num2str(k)));
    p.fileax(k)=ax;
    
    [imstack, p.roi{k}, p.pixelsize{k},settings3D]=readbeadimages(filelist{k},p);
    if ~isempty(settings3D)
        p.settings_3D=settings3D;
    end
    
%     if isfield(p,'smap') && p.smap %at some point we have to overcome this!
%         try
%              r=imageloaderAll(filelist{k},[],p.smappos.P);
% 
%              imstack=r.getmanyimages(1:r.metadata.numberOfFrames,'mat');
%              p.roi{k}=r.metadata.roi;
%              p.pixelsize{k}=r.metadata.cam_pixelsize_um;
%              r.close;
%         catch err
%             err
% %             imstack=readfile_ome(filelist{k});
% %             p.roi{k}=p.smappos.roi{k};
% %             p.pixelsize{k}=p.smappos.pixelsize{k};
%             imstack=readfile_tif(filelist{k});
%             p.roi{k}=[0 0 size(imstack,1) size(imstack,2)]; %check x,y
%         end
%         if isempty(imstack)
%             disp('using simple reader')
%             imstack=readfile_tif(filelist{k});
%         end
%           
%     else
%         imstack=readfile_tif(filelist{k});
%         p.roi{k}=[0 0 size(imstack,1) size(imstack,2)]; %check x,y
%     end
    
    if isfield(p,'settings_3D') %calibration file: cut out and mirror already here!
        imstack=cutoutchannels(imstack,p.settings_3D);
        p.roi{k}=[0 0 size(imstack,1) size(imstack,2)]; %check x,y
    end
    if p.emgain
        imstack=imstack(:,end:-1:1,:);
    end
       
    if isfield(p,'framerangeuse')
        imstack=imstack(:,:,p.framerangeuse(1):p.framerangeuse(end));
    end
    
    imstack=imstack-min(imstack(:)); %fast fix for offset;
    
%     imageslicer(imstack)%%%%%%%%XXXXXXX
    
    mim=max(imstack,[],3);
    mim=mean(imstack,3);
    mim=filter2(h,mim);

    imagesc(ax,mim);
    axis(ax,'equal');
    axis(ax,'off')
    title(ax,'Maximum intensity projection')
    if isfield(p,'beadpos') %passed on externally
        maxima=round(p.beadpos{k});
    else
        if isfield(p,'roimask')&&~isempty(p.roimask)
            mim=mim.*(p.roimask);
        end
        maxima=maximumfindcall(mim);
        %XXXXXXXXXXX
%         maxima(:,1)=(maxima(:,1)+2*round(rand(size(maxima,1),1)))-1; % for testing if positions match
        int=maxima(:,3);
        try
        mimc=mim(roisize:end-roisize,roisize:end-roisize);
        mmed=myquantile(mimc(:),0.3);
        imt=mimc(mimc<mmed);
            sm=sort(int);
        mv=mean(sm(end-5:end));
%         cutoff=mean(imt(:))+max(2.5*std(imt(:)),(mv-mean(imt(:)))/10);
        cutoff=mean(imt(:))+max(2.5*std(imt(:)),(mv-mean(imt(:)))/15);
        catch
            cutoff=myquantile(mimc(:),.95);
        end
        if isnan(cutoff)
            cutoff=myquantile(mimc(:),.95);
        end
        cutoff=cutoff*p.cutoffrel;
        if any(int>cutoff)
            maxima=maxima(int>cutoff,:);
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
        indgoodr=maxima(:,1)>p.xrange(1)&maxima(:,1)<p.xrange(end)&maxima(:,2)>p.yrange(1)&maxima(:,2)<p.yrange(end);
        if isfield(p,'settings_3D')
           w=p.settings_3D.width4pi;
           xm=mod(maxima(:,1),w);
           xm(xm>w/2)=xm(xm>w/2)-w;
           indgoodb(abs(xm)<p.mindistance/2)=false;
       end
       h=size(imstack,1);
       ym=maxima(:,2);
       ym(ym>h/2)=ym(ym>h/2)-h;
       indgoodb(abs(ym)<p.mindistance/2)=false;
%         maxima=maxima(indgoodb,:);
         maxima=maxima(indgoodb&indgoodr,:);
    end 
    
  
    if p.isglobalfit
        transform=p.transformation;
        %calculate in nm on chip (reference for transformation)
        maximapx=(maxima(:,1:2)+p.roi{k}([1 2]));
        
%         if isfield(transform.tinfo,'units') &&strcmp(transform.tinfo.units,'pixels')
%             pixfac=[1 1];
%         else
%             pixfac=[p.pixelsize{k}(1)*1000 p.pixelsize{k}(end)*1000];
%         end
%             maximanm(:,1)=maximanm(:,1)*pixfac(1);
%             maximanm(:,2)=maximanm(:,2)*pixfac(end);

        %transform reference to target

            indref=transform.getPart(1,maximapx);
            maximaref=maxima(indref,1:2);
            maximaall=transform.transformToTargetAll(maximaref);
%             [x,y]=transform.transformCoordinatesFwd(maximapx(indref,1),maximapx(indref,2));
        %     [x,y]=transform.transformCoordinatesFwd(maximanm(indref,2),maximanm(indref,1));

%             maximatargetf=[];
            maximaall(:,1,:)=maximaall(:,1,:)-p.roi{k}(1);
            maximaall(:,2,:)=maximaall(:,2,:)-p.roi{k}(2);
%             maximatargetf(:,1)=x/pixfac(1)-p.roi{k}(1);
%             maximatargetf(:,2)=y/pixfac(end)-p.roi{k}(2);
    %     maximatargetf(:,2)=x/p.smappos.pixelsize{k}(1)/1000-p.smappos.roi{k}(1);
    %     maximatargetf(:,1)=y/p.smappos.pixelsize{k}(end)/1000-p.smappos.roi{k}(2);


%         if 1 %for testing
%             maximatargetf(:,1)=maximatargetf(:,1)+1;
%             maximatargetf(:,2)=maximatargetf(:,2)+0.5;
%         maximatargetfm(:,1)=maximatargetf(:,1)-0.1+2;
%         maximatargetfm(:,2)=maximatargetf(:,2)+0.1;
%         maximatar=round(maximatargetfm); 
%         dxy=maximatargetf-maximatar;  
%         else 
            maximaallr=round(maximaall); 
%             maximaallr(:,1,3)=maximaallr(:,1,3)+2;
            
           
            dxy=maximaall-maximaallr;       

%         end
     
    else
        
        maximaref=maxima;
        maximatar=maxima;
        dxy=zeros(size(maximatar));
        maximaallr=maxima;
    end
   
    numframes=size(imstack,3);
    bind=length(b)+size(maximaref,1);

    for l=1:size(maximaref,1)
        b(bind).loc.frames=(1:numframes)';
        b(bind).loc.filenumber=zeros(numframes,1)+k;
        b(bind).filenumber=k;
        b(bind).posall=squeeze(maximaallr(l,:,:));
        b(bind).pos=squeeze(maximaallr(l,:,1));
%         b(bind).postar=maximatar(l,1:2);
        b(bind).shiftxy=dxy(l,:,:);
        try
            if p.isglobalfit %many channels present
                for ch=1:size(maximaallr,3)
                    b(bind).stack.image(:,:,:,ch)=imstack(b(bind).posall(2,ch)+rsr,b(bind).posall(1,ch)+rsr,:);
%                 b(bind).stack.imagetar=imstack(b(bind).postar(2)+rsr,b(bind).postar(1)+rsr,:);
                end
            else
                b(bind).stack.image(:,:,:)=imstack(b(bind).posall(2)+rsr,b(bind).posall(1)+rsr,:);
            end
            b(bind).stack.framerange=1:numframes;
            b(bind).isstack=true;
            
        catch err
            b(bind).isstack=false;
%             err
        end
        
        b(bind).roi=p.roi{k};
        bind=bind-1;
    end
    fmax=max(fmax,numframes);
     hold (ax,'on')
    if p.isglobalfit
        for ch=1:size(maximaall,3)
            plot(ax,maximaall(:,1,ch),maximaall(:,2,ch),'ko',maximafound(:,1),maximafound(:,2),'r.') 
        end
    else
        plot(ax,maxima(:,1),maxima(:,2),'ko',maximafound(:,1),maximafound(:,2),'r.')
    end
    hold (ax,'off')
    drawnow
end
indgoodbead=[b(:).isstack];
b=b(indgoodbead);

% plot
p.fminmax=[1 fmax];

p.pathhere=fileparts(filelist{1});
end


