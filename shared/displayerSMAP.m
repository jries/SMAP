function [imout,sr_imagehandle]=displayerSMAP(layers,p)
% Combines rendered images from different channels and adds scale and color bars.
if nargin==0
    %input parameters
    imout={'sr_plotcomposite','sr_layerson','sr_axes','sr_sizeRecPix','roihandle','sr_pixrec','rotationangle','sr_pos','sr_size','sr_layersseparate','layernames','sr_plotlayernames','sr_plotscalebar','sr_colorbarthickness','sr_lutwhite'};
    return          
end

layersnext=isfield(p,'sr_layersseparate')&&~isempty(p.sr_layersseparate)&&p.sr_layersseparate;
        plotcomposite=isfield(p,'sr_plotcomposite')&&~isempty(p.sr_plotcomposite)&&p.sr_plotcomposite; % also plot composite

tiffthere=0;
txtN='';
allnext=[];



for k=1:length(layers)

    if p.sr_layerson(k)
        if isfield(layers(k),'images') && ~isempty(layers(k).images)
%             if show==0
%                 si=size(layers(k).images.finalImages.image);
%                  imall=zeros(si(1)+rim,si(2)+rim,3);
%                  mask=zeros(si(1)+rim,si(2)+rim,3);
%                  imtiff=zeros(si(1)+rim,si(2)+rim,3);
% %                  imref=layers(k).images.finalImages;
% 
%                 show=1;
%             end
         fi=layers(k).images.finalImages;
         if isempty(fi.image)
             continue
         end
         if fi.istiff==1
             if exist('imtiff','var')
                imtiff=imtiff+fi.image;
             else
                 imtiff=fi.image;
             end
                tiffthere=1;
         else
             if exist('imall','var')
                 if ~all(size(imall)==size(fi.image))%s&&~layersnext
                    simh=size(imall); 
                    fi.image=imresize(fi.image,simh(1:2));
                    fi.mask=imresize(fi.mask,simh(1:2));
                 end
                 imall=imall+fi.image;
                 mask=mask+fi.mask;
             else
                 imall=fi.image;
                 mask=fi.mask;
             end
              txtN=[txtN 'N'  num2str(k) '=' shortnumber(fi.numberOfLocs) ', '];
         end
         if layersnext
             if fi.istiff && exist('imall','var')
                 if ~all(size(imall)==size(fi.image))%s&&~layersnext
                    simh=size(imall); 
                    fi.image=imresize(fi.image,simh(1:2));
                 end
             end
             s=size(fi.image);
             if s(2)>s(1)*1.3 
                 vertnext=true;
                 fi.image(end,:,:)=.5;
                 allnext=vertcat(fi.image,allnext);
             else
                 fi.image(:,end,:)=.5;
                 allnext=horzcat(allnext,fi.image);
                 vertnext=false;
             end
         end

        end
    end
end


%make color bars
if ~exist('imall','var')&&~tiffthere
    imout=[];
    sr_imagehandle=[];
    return
end

if ~exist('imall','var')
    imfinal=imtiff;
elseif tiffthere
    sm=size(mask);
    if min(sm(1:2))>4
        mask(1:4,:,:)=1;
        mask(:,1:4,:)=1;
        mask(end-4:end,:,:)=1;
        mask(:,end-4:end,:)=1;
    end
    mask(mask>1)=1;
    if any(size(imtiff) ~=size(imall))
        sima=size(imall);
        imtiff=imresize(imtiff,sima(1:2));
    end
    
    imfinal=mask.*imall+(1-mask).*imtiff;
else
    imfinal=imall;
end

compimage=imfinal;    
%rotate
if isfield(p,'rotationangle')&&~isempty(p.rotationangle)&&p.rotationangle~=0
    imfinal=imrotate(imfinal,p.rotationangle,'nearest','crop');
end
if ~isfield(p,'sr_colorbarthickness')
    p.sr_colorbarthickness=4;
end
for k=1:(length(layers))
    if p.sr_layerson(k)&&~isempty(layers(k).images)
%         if k<=4 && (isempty(p.sr_colorbarthickness) ||  p.sr_colorbarthickness>0)
%         imfinal=addcolorbar(imfinal,layers(k).images.finalImages.lut,k,p.sr_colorbarthickness);
%         end
        rangexplot=layers(k).images.finalImages.rangex;
         rangeyplot=layers(k).images.finalImages.rangey;
    end
end


if layersnext
    if plotcomposite
         if vertnext
             allnext=vertcat(imfinal,allnext);
         else
             allnext=horzcat(allnext,imfinal);
         end
    end
    
     imfinal=allnext;
     nlayer=sum(p.sr_layerson)-1+plotcomposite;
     if vertnext
         rangeyplot(2)=rangeyplot(2)+nlayer*(rangeyplot(2)-rangeyplot(1));
     else
        rangexplot(2)=rangexplot(2)+nlayer*(rangexplot(2)-rangexplot(1));
     end
end

if (isfield(p,'addscalebar') && ~ p.addscalebar) || (isfield(p,'sr_plotscalebar') && ~isempty(p.sr_plotscalebar) && ~ p.sr_plotscalebar) 
     lennm=0;
else
    [imfinal,lennm]=addscalebar(imfinal,p.sr_pixrec(1));
   
end
    
colorpos=1;
for k=1:(length(layers))
    if p.sr_layerson(k)&&~isempty(layers(k).images)
        if k<=4 && (isempty(p.sr_colorbarthickness) ||  p.sr_colorbarthickness>0)
        imfinal=addcolorbar(imfinal,layers(k).images.finalImages.lut,colorpos,p.sr_colorbarthickness);
        colorpos=colorpos+1;
        end
%         rangexplot=layers(k).images.finalImages.rangex;
%          rangeyplot=layers(k).images.finalImages.rangey;
    end
end
    
if isfield(p,'sr_lutwhite') && ~isempty(p.sr_lutwhite) && p.sr_lutwhite
%     imfinal=invertwhite(imfinal,0.0);
    imfinal=1-imfinal;
end

    if isfield(p,'sr_axes')&&~isempty(p.sr_axes)&&ishandle(p.sr_axes)&&~isempty(rangexplot)&&~isempty(rangeyplot)
        sr_imagehandle=image(rangexplot/1000,rangeyplot/1000,imfinal,'Parent',p.sr_axes,'Pickable','none','HitTest','off');
%                     plotovim=1;
        title(p.sr_axes,txtN)
        set(p.sr_axes,'Xlim',rangexplot/1000)
        set(p.sr_axes,'Ylim',rangeyplot/1000)
        set(p.sr_axes,'YDir','reverse')
        axis(p.sr_axes,'equal')
        xlabel(p.sr_axes,'x (µm)')
        ylabel(p.sr_axes,'y (µm)')
        p.sr_axes.HitTest='on';
        p.sr_axes.PickableParts='all';
%         axes(p.sr_axes) %bring to forground
        imout.handle=sr_imagehandle;
        drawnow limitrate nocallbacks
        
        dxy=0;
        fontsize=15;
        extent=0;
        if isfield(p,'sr_plotlayernames')&&~isempty(p.sr_plotlayernames)&&p.sr_plotlayernames
             for k=1:(length(layers))
                 if p.sr_layerson(k)&&~isempty(layers(k).images)

        
                     lut=layers(k).images.finalImages.lut;

                     
                     if layersnext
                        dx=(layers(k).images.srimage.rangex(2)-layers(k).images.srimage.rangex(1))/1000*dxy;
                        dy=0;
                     else
                        dx=0;
                        dy=dxy;
                     end
                     px=double(layers(k).images.srimage.rangex(1)/1000+dx+p.sr_pixrec/1000*5);
                     py=double(layers(k).images.srimage.rangey(1)/1000);
                     lutm=mean(lut,1);lutm=lutm/max(lutm);
                     th=text(p.sr_axes,px,py,p.layernames{k},'Color',lutm,'FontSize',fontsize,'BackgroundColor','k','Units','data');
%                      th.Units='pixels';
                     th.Position(2)=th.Position(2)+th.Extent(4)*(dy+.6)*1.2;
                     
                     dxy=dxy+1;
                 end
             end
             if layersnext && plotcomposite
                 dx=(layers(k).images.srimage.rangex(2)-layers(k).images.srimage.rangex(1))/1000*dxy;
                 px=layers(k).images.srimage.rangex(1)/1000+dx+p.sr_pixrec/1000*5;
                 py=layers(k).images.srimage.rangey(1)/1000;
                 th=text(p.sr_axes,px,py,'composite','Color','w','FontSize',fontsize,'BackgroundColor','k','Units','data');
                 th.Position(2)=th.Position(2)+th.Extent(4)*(dy+.6)*1.2;
             end
            
        end
    else
        sr_imagehandle=[];
    end
    

    
    imout.image=imfinal;
    imout.composite=compimage;
    imout.rangex=rangexplot/1000;
    imout.rangey=rangeyplot/1000;
    imout.scalebarnm=lennm;

end


function imout=addcolorbar(imin,lut,layer,thickness)
if isempty(imin)
    imout=imin;
    return
end
s=size(imin);
if layer==1||layer==3
    l=s(1);
else
    l=s(2);
end
if layer<4
    lld=0;
else
    lld=2;
end

ll=thickness+lld;
rim=10;

x=ceil((1:l-rim)/(l-rim)*255);
y=1:ll;
[X,~]=meshgrid(x,y);
lutim=ind2rgb(uint8(X),lut);
imout=imin;
switch layer
    case 3
        imout(rim/2+1:end-rim/2,1:ll,:)=permute(lutim,[2 1 3]);
    case 2
%         imout(end-3:end,rim/2+1:end-rim/2,:)=lutim;
        imout(1:ll,rim/2+1:end-rim/2,:)=lutim;
    case 1
        imout(rim/2+1:end-rim/2,end-ll+1:end,:)=permute(lutim(:,end:-1:1,:),[2 1 3]);
    case 4
%         imout(1:6,rim/2+1:end-rim/2,:)=lutim;
        imout(end-ll+1:end,rim/2+1:end-rim/2,:)=lutim;
end

end

function [imin,lennm]=addscalebar(imin,pixrec,fac)
sim=size(imin);
if nargin<3
    fac=1;
end
   lennm=10^floor(log10(sim(2)*pixrec*0.7*fac));
   lenpix=round(lennm/pixrec);
   thickness=max(2,round(sim(2)/200));
   if lenpix<sim(2)-12&&sim(1)>4
   imin(end-thickness-2:end-1,end-10-lenpix:end-9,:)=0;
   imin(end-thickness-1:end-2,end-10-lenpix+1:end-10,:)=1;
   end
end
