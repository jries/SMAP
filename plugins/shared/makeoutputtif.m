function outimage=makeoutputtif(obj,p)
if nargin==0  %default, can be used to initialize gui control
    outimage={'rendered image with scalebar','rendered image','layer 1-3 as RGB','layer1 as grayscale'};
    return
end
if ~isfield(p,'addscalebar')
    p.addscalebar=false;
end
 switch p.outputformat.selection
    case 'rendered image with scalebar'
        srimage=obj.getPar('sr_image');
        outimage=uint8(srimage.image*255);
        p.addscalebar=false;
%                        mij.createColor(title,outimage,true);
    case 'rendered image'
        srimage=obj.getPar('sr_image');
        outimage=uint8(srimage.composite*255);
        maxv=255;
%                        mij.createColor(title,outimage,true);
    case 'layer 1-3 as RGB'
        s1=obj.locData.layer(1).images.srimage.image;
        sizes=size(s1);
        if length(sizes)>2
            s1=sum(s1,3);
        end
        outimage=zeros(sizes(1),sizes(2),3);
        outimage(:,:,1)=s1/obj.locData.layer(1).images.finalImages.imax;
        if length(obj.locData.layer)>1
             s2=obj.locData.layer(2).images.srimage.image;
              if length(size(s2))>2
                s2=sum(s2,3);
              end
              if size(s2)==size(s1)
                  outimage(:,:,2)=s2/obj.locData.layer(2).images.finalImages.imax;
              end
        end
        if length(obj.locData.layer)>2
             s3=obj.locData.layer(3).images.srimage.image;
              if length(size(s3))>2
                s3=sum(s3,3);
              end
              if size(s3)==size(s1)
                  outimage(:,:,3)=s3/obj.locData.layer(3).images.finalImages.imax;
              end
        end
        outimage(outimage>1)=1;
        outimage=uint8(outimage*255);
%         outimage=uint8(outimage/max(outimage(:))*255);
        maxv=255;
%                        mij.createColor(title,outimage,true);

    case 'layer1 as grayscale'
        s1=obj.locData.layer(1).images.srimage.image;
        sizes=size(s1);
        if length(sizes)>2
            s1=sum(s1,3);
        end
        outimage=s1;
        outimage=uint16(outimage/max(outimage(:))*(2^16));
        maxv=2^16;
end
%                        mij.createImage(title,outimage,true);
if p.addscalebar
    srimage=obj.getPar('sr_image');
    pixelsize=(srimage.rangey(2)-srimage.rangey(1))/size(srimage.composite,1)*1000;
    scalebarlength=round(srimage.scalebarnm/pixelsize);
    outimage(end-4:end,end-scalebarlength-10+1-1:end-10+1,:)=0;
    outimage(end-3:end-1,end-scalebarlength-10+1:end-10,:)=maxv;
end
end

 