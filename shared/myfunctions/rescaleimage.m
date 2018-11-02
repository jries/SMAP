function [imout,maxi]=rescaleimage(imin,gamma,quant)

s=size(imin);

%  imout=im.^gamma;

switch length(s)
    case 2 
        intim=(imin);
         intim=intim/max(intim(:));
         cal=intim;
         
         imout=imin./cal.*intim.^gamma;

         imout(isnan(imout))=0;


          mini=min(imout(:));
        maxi=myquantilefast(imout(:),quant,30/(1-quant));
        if maxi==mini
            maxi=mini+1
        end

        imin=(imout-mini)/(maxi-mini)*1; %rescale

        cutoff=1.3;
        intmax=sum(imin,3);

            iminc=imin;
            iminc(intmax>cutoff)=iminc(intmax>cutoff)./intmax(intmax>cutoff)*cutoff;
          imout=iminc;


         imout(imout>1)=1;
                
        
    case 3 %normal image XXXXXXXXXXwhat if BW movie???
         intim=sum(imin,3);
         intim=intim/max(intim(:));
         cal=intim;
         imout=0*imin;
         for k=1:3
            imout(:,:,k)=imin(:,:,k)./cal.*intim.^gamma;
         end
         imout(isnan(imout))=0;
         

         for k=1:3  
            
            imhere=imout(:,:,k); 
            
            mini=min(imhere(:));
            imhere=imhere-mini;
            
            imquant=imhere(end/2-end/4:end/2+end/4,end/2-end/4:end/2+end/4);
            maxi=myquantilefast(imquant(:),quant,30/(1-quant));
            imout(:,:,k)=(imhere/maxi);
         end
         
%what is this???
%         cutoff=1.3;
%         intmax=sum(imin,3);
%          for k=1:3
%             iminc=imin(:,:,k);
%             iminc(intmax>cutoff)=iminc(intmax>cutoff)./intmax(intmax>cutoff)*cutoff;
%           imout(:,:,k)=iminc;
%          end

         imout(imout>1)=1;
 
    case 4 %movie
 
         intim=sum(imin,3);
         intim=intim/max(intim(:));
         cal=intim;
         imout=0*imin;
         for k=1:3
         imout(:,:,k,:)=imin(:,:,k,:)./cal.*intim.^gamma;
         end
         imout(isnan(imout))=0;


          mini=min(imout(:));
        maxi=myquantilefast(imout(:),quant,30/(1-quant));
                if maxi==mini
                    maxi=mini+1
                end
        imin=(imout-mini)/(maxi-mini)*1; %rescale

        cutoff=1.3;
        intmax=sum(imin,3);
         for k=1:3
            iminc=imin(:,:,k,:);
            iminc(intmax>cutoff)=iminc(intmax>cutoff)./intmax(intmax>cutoff)*cutoff;
          imout(:,:,k,:)=iminc;
         end

         imout(imout>1)=1;
end
 
