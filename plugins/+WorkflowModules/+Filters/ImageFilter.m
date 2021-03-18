classdef ImageFilter<interfaces.WorkflowModule
%     Filters the image before peak finding. Included are a Gaussian filter
%     (recommended for beads and if the background is determined
%     independentely), a Difference-of-Gaussian filter  (works well on most
%     SMLM data) and filtering with an experimental PSF (recommended for
%     non-compact PSFs, e.g. the double-helix PSF.
    properties
        filterkernel
        preview
        filterkernelPSF
        offset=0;
        scmosvariance
        scmosfiltered
    end
    methods
        function obj=ImageFilter(varargin)
            obj@interfaces.WorkflowModule(varargin{:});
            obj.outputParameters={'loc_loc_filter_sigma'};
            
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function prerun(obj,p)
            obj.scmosvariance=[];
            p=obj.getAllParameters;
            fs=p.loc_loc_filter_sigma(1);
            if length(p.loc_loc_filter_sigma)>1
                obj.offset=p.loc_loc_filter_sigma(2);
            end
            if fs>0
                h=fspecial('gaussian', max(5,ceil(3.5/2*fs)*2+1), fs);
            else
                h=1;
            end
            switch p.filtermode.Value
                case 2 %GAuss
                    obj.filterkernel=h;

                case 3
                        PSF=obj.filterkernelPSF.PSF;
                        fig=nanmean(PSF,3);
                        rs=obj.getPar('loc_ROIsize');
                        rss=round((rs-1)/2);
                        midp=round(size(fig)/2);
                        figs=fig(midp(1)-rss:midp(1)+rss,midp(2)-rss:midp(2)+rss);
                        kernelPSF=figs;%-nanmin(figs(:));
                        kernelPSF=kernelPSF/nansum(kernelPSF(:));
                     
                    if isempty(obj.filterkernelPSF)
                        error('for PSF filtering you need to load a PSF model _3Dcal.mat before fitting')
                    end
                    if fs>0
                        obj.filterkernel=filter2(h,kernelPSF);
                    else
                        obj.filterkernel=kernelPSF;
                    end
                case 4 %exp PSF, MIP
                        rs=obj.getPar('loc_ROIsize');
                        rss=round((rs-1)/2);
                        
                        PSF=obj.filterkernelPSF.PSF;
                        midp=round(size(PSF)/2);
                        dz=obj.filterkernelPSF.dz;
                        if numel(p.zrange)>1
                            z=round(p.zrange/dz+size(PSF,3)/2);
                            z(z<1)=1;z(z>size(PSF,3))=size(PSF,3);
                        else
                            z=round(linspace(1,size(PSF,3),max(2,p.zrange)));
                        end
                        psfstack=zeros(rs,rs,length(z)-1);
                        normp=PSF(midp(1)-rss:midp(1)+rss,midp(2)-rss:midp(2)+rss,round(size(PSF,3)/2));
                        normf=nansum(normp(:));
%                         minph=nanmin(normp(:));
                        for k=1:length(z)-1
                            ph=nanmean(PSF(midp(1)-rss:midp(1)+rss,midp(2)-rss:midp(2)+rss,z(k):z(k+1)-1),3);
%                             ph=ph-min(ph(:));
                            psfstack(:,:,k)=ph/normf;
%                             psfstack(:,:,k)=ph/nansum(ph(:));
                        end
                        obj.filterkernelPSF.fitpsf=psfstack;
                        if fs>0
                            obj.filterkernel=h;
                        else 
                            obj.filterkernel=[];
                        end
                        
                case 1 %DoG
                    rsize=max(ceil(6*fs-1),3);
                    if p.correctsCMOS
                        obj.filterkernel=[];
                        obj.filterkernel(:,:,1)=fspecial('gaussian',rsize,fs);
                        obj.filterkernel(:,:,2)=fspecial('gaussian',rsize,max(1,2.5*fs));
                    else
                        hdog=fspecial('gaussian',rsize,fs)-fspecial('gaussian',rsize,max(1,2.5*fs));
                        obj.filterkernel=hdog;
                    end
            end
            obj.preview=obj.getPar('loc_preview');
        end
        function dato=run(obj,data,p)
            if isempty(data.data)
                dato=data;
                return
            end
            dato=data;%.copy;
            if p.filtermode.Value==1&&~isempty(data.data)&&obj.offset==0
                offset=min(data.data(:,1));
            else 
                offset=obj.offset;
            end
            
            if p.correctsCMOS && isempty(obj.scmosvariance)  %initialize scmose variance map
                vmap=obj.getPar('cam_varmap');
                vsmallest=quantile(vmap(:),.25); %dont overweigh pixels with too low variance (broken pixels)
                vmap(vmap<vsmallest)=vsmallest;
%                 vlargest=100*median(vmap(:));
%                 vmap(vmap>vlargest)=vlargest;
                obj.scmosvariance=vmap;
                switch p.filtermode.Value
                    case 1
                        obj.scmosfiltered(:,:,1)=filter2(obj.filterkernel(:,:,1),1./obj.scmosvariance);
                        obj.scmosfiltered(:,:,2)=filter2(obj.filterkernel(:,:,2),1./obj.scmosvariance);
                    case 2
                        obj.scmosfiltered=filter2(obj.filterkernel,1./obj.scmosvariance);
                    case 3
                        obj.scmosfiltered=conv2(1./obj.scmosvariance,obj.filterkernel,'same');
                    case 4
                        disp('scmose noise correction in peak finder not implemented for MIP');
                end
            end
            
            if p.correctsCMOS && p.filtermode.Value<4
               switch p.filtermode.Value
                   case 1
                       imf=filter2(obj.filterkernel(:,:,1),(data.data-offset)./obj.scmosvariance)./obj.scmosfiltered(:,:,1)-filter2(obj.filterkernel(:,:,2),(data.data-offset)./obj.scmosvariance)./obj.scmosfiltered(:,:,2);
                    case 2
                        imf=filter2(obj.filterkernel,(data.data-offset)./obj.scmosvariance)./obj.scmosfiltered;
                    case 3
                        dat=(data.data/2).^2-0.375;
                        imf=conv2(dat./obj.scmosvariance,obj.filterkernel,'same').*obj.scmosfiltered;     
                end
            else
                switch p.filtermode.Value
                    case {1,2}
                        imf=filter2(obj.filterkernel,data.data-offset);
                    case 3
                        dat=(data.data/2).^2-0.375;
                        imf=conv2(dat,obj.filterkernel,'same');
                    case 4
                        h=obj.filterkernel;
                        psfstack=obj.filterkernelPSF.fitpsf;
                        imf=-inf;
                        dat=data.data;
                        for k=1:size(psfstack,3)
                            dat=(data.data/2).^2-0.375;
                            imh=conv2(dat,psfstack(:,:,k),'same');
                            imf=max(imh,imf);
                        end
                        if ~isempty(obj.filterkernel)
                            imf=filter2(obj.filterkernel,imf);
                        end
                        imf=(2*sqrt(imf+0.3750));     
                end
            end
            
            if obj.preview
                drawimage(obj,data.data-offset,imf)   
            end
           

            dato.data=(imf);
        end
    end
end


function drawimage(obj,imnorm,imf)
if isempty(imf)
    return
end

img=(imnorm/2).^2-0.375;
imgbg=(imf/2).^2-0.375;

imbg=((imnorm-imf)/2).^2-0.375;

% outputfig=obj.getPar('loc_outputfig');
% if ~isvalid(outputfig)
%     outputfig=figure(209);
%     obj.setPar('loc_outputfig',outputfig);
% end

% outputfig.Visible='on';
% draw=~isempty(imnorm);
% switch obj.getPar('loc_previewmode').Value
%     case 1 %image-bg
%         imd=imbg;
%     case 2%image
%         imd=img;
%     case 3 %norm
%         imd=imf;
%     case 4 %bg
%         imd=imgbg;
%     otherwise 
%         draw=false;
% end
        
% if draw
% figure(outputfig)
% hold off
% imagesc(imd);
% colormap jet
% colorbar;
% axis equal
% hold on
% end
obj.setPar('preview_filtered',imf);
obj.setPar('preview_background',imgbg);
obj.setPar('preview_image_background',imbg);
end


function loadPSF_callback(object,b,obj)

p=(obj.getPar('lastSMLFile'));
if isempty(p)
    p=obj.getPar('loc_fileinfo');
    p=p.basefile;
end
if ~isempty(p)
    p=fileparts(p);
end
[f,p]=uigetfile([p filesep '*.mat']);

if f  
    l=load([p f]);
    if isfield(l,'SXY')
        l=l.SXY(1);
    end
    if isfield(l,'PSF')
        PSF=l.PSF{1};
        dz=l.cspline.dz;
    else
        disp('PSF not found')
    end
    PSF=PSF-nanmin(PSF(:));
%     PSF=(2*sqrt(PSF-nanmin(PSF(:))+0.3750)); % image normalize
    obj.filterkernelPSF.PSF=PSF;
    obj.filterkernelPSF.dz=dz;
    
end
end

function pard=guidef(obj)
p(1).value=1;p(1).on={};p(1).off={'loadPSF','text2','zrange'};
p(2)=p(1);p(2).value=2;
p(3).value=3;p(3).on={'loadPSF'};p(3).off={'text2','zrange'};
p(4).value=4;p(4).on={'loadPSF','text2','zrange'};p(4).off={};

pard.filtermode.object=struct('Style','popupmenu','String',{{'DoG','Gauss: ','mean PSF','MIP PSF'}},'Callback',{{@obj.switchvisible,p}});
pard.filtermode.position=[1,1];
pard.filtermode.Width=0.8;
pard.filtermode.Optional=true;

pard.text.object=struct('Style','text','String','s:');
pard.text.position=[1,1.8];
pard.text.Width=0.15;
pard.text.Optional=true;

pard.loc_loc_filter_sigma.object=struct('Style','edit','String','1.2','Visible','on');
pard.loc_loc_filter_sigma.position=[1,1.95];
pard.loc_loc_filter_sigma.Width=.3;
pard.loc_loc_filter_sigma.TooltipString=sprintf('Sigma (in camera pixels) for a Gaussian filter which is applied after background correction and before peak finding. \n Typical size of PSF in pixels, eg 1 (range: 0.5-5) ');

pard.loadPSF.object=struct('Style','pushbutton','String','load','Callback',{{@loadPSF_callback,obj}},'Visible','off');
pard.loadPSF.position=[1,2.2];
pard.loadPSF.Width=.4;

pard.text2.object=struct('Style','text','String','z:','Visible','off');
pard.text2.position=[1,2.6];
pard.text2.Width=0.2;
pard.text2.Optional=true;

pard.zrange.object=struct('Style','edit','String','5');
pard.zrange.position=[1,2.8];
pard.zrange.Width=1;
pard.zrange.TooltipString=sprintf('For using experimetnal PSF for peak finding: at which z-positions to probe the PSF for correlation');


pard.correctsCMOS.object=struct('Style','checkbox','String','correct sCMOS variance');
pard.correctsCMOS.position=[2,1.25];
pard.correctsCMOS.Width=1.75;
pard.correctsCMOS.TooltipString=sprintf('Correct sCMOS variance');

pard.plugininfo.type='WorkflowModule';
pard.loc_loc_filter_sigma.Optional=true;
pard.plugininfo.description='Filters the image before peak finding. Included are a Gaussian filter (recommended for beads and if the background is determined independentely), a Difference-of-Gaussian filter  (works well on most SMLM data) and filtering with an experimental PSF (recommended for non-compact PSFs, e.g. the double-helix PSF.';
pard.text.TooltipString=pard.loc_loc_filter_sigma.TooltipString;
end