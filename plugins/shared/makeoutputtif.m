function outimage=makeoutputtif(obj,p)
if nargin==0  %default, can be used to initialize gui control
    outimage={'rendered image with scalebar','rendered image','layer 1-3 as RGB','layer1 as grayscale','layers plus composite'};
    return
end
if ~isfield(p,'addscalebar')
    p.addscalebar=false;
end
 switch p.outputformat.selection
    case 'rendered image with scalebar'
        srimage=obj.getPar('sr_image');
        outimage{1}=uint8(srimage.image*255);
        p.addscalebar=false;
%                        mij.createColor(title,outimage,true);
    case 'rendered image'
        srimage=obj.getPar('sr_image');
        outimage{1}=uint8(srimage.composite*255);
        maxv=255;
%                        mij.createColor(title,outimage,true);
    case 'layer 1-3 as RGB'
        s1=obj.locData.layer(1).images.srimage.image;
        sizes=size(s1);
        if length(sizes)>2
            s1=sum(s1,3);
        end
        outimage{1}=zeros(sizes(1),sizes(2),3);
        outimage{1}(:,:,1)=s1/obj.locData.layer(1).images.finalImages.imax;
        if length(obj.locData.layer)>1
             s2=obj.locData.layer(2).images.srimage.image;
              if length(size(s2))>2
                s2=sum(s2,3);
              end
              if size(s2)==size(s1)
                  outimage{1}(:,:,2)=s2/obj.locData.layer(2).images.finalImages.imax;
              end
        end
        if length(obj.locData.layer)>2
             s3=obj.locData.layer(3).images.srimage.image;
              if length(size(s3))>2
                s3=sum(s3,3);
              end
              if size(s3)==size(s1)
                  outimage{1}(:,:,3)=s3/obj.locData.layer(3).images.finalImages.imax;
              end
        end
        outimage{1}(outimage{1}>1)=1;
        outimage{1}=uint8(outimage*255);
%         outimage=uint8(outimage/max(outimage(:))*255);
        maxv=255;
%                        mij.createColor(title,outimage,true);

    case 'layer1 as grayscale'
        s1=obj.locData.layer(1).images.srimage.image;
        sizes=size(s1);
        if length(sizes)>2
            s1=sum(s1,3);
        end
        outimage{1}=s1;
        outimage{1}=uint16(outimage/max(outimage(:))*(2^16));
        maxv=2^16;
     case 'layers plus composite'
         layers=find(obj.getPar('sr_layerson'));        
         for k=1:length(layers)
             sm=displayerSMAP(obj.locData.layer(layers(k)),p);
             outimage{k}=uint8(sm.composite*255);
         end
         srimage=obj.getPar('sr_image');
         outimage{k+1}=uint8(srimage.composite*255);
         maxv=255;
end
%                        mij.createImage(title,outimage,true);
if p.addscalebar
    srimage=obj.getPar('sr_image');
    pixelsize=(srimage.rangey(2)-srimage.rangey(1))/size(srimage.composite,1)*1000;
    scalebarlength=round(srimage.scalebarnm/pixelsize);
    for k=1:length(outimage)
        outimage{k}(end-4:end,end-scalebarlength-10+1-1:end-10+1,:)=0;
        outimage{k}(end-3:end-1,end-scalebarlength-10+1:end-10,:)=maxv;
    end
end
end

 