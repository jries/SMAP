function [imphot,varmap]=correctcameraimages(imstack,metadata,correctemmirror,correctscmos)
if nargin<3
    correctemmirror=true;
end
if nargin<4
    correctscmos=true;
end
imscmos=[];
imgp=single(makepositive(imstack));

if ~metadata.EMon
    emgain=1;
else 
    emgain=metadata.emgain;
end
adu2phot=(metadata.conversion/emgain);
            
            
if correctscmos %apply offset and brightfield correction
   roi=metadata.roi;
   [gainmap,offsetmap,varmap]=makegainoffsetCMOS(metadata.correctionfile,metadata.exposure);
   if ~isempty(gainmap)
       imscmos=(single(imgp)-offsetmap(roi(1)+1:roi(1)+roi(3),roi(2)+1:roi(2)+roi(4)))... %*obj.adu2phot...
           .*gainmap(roi(1)+1:roi(1)+roi(3),roi(2)+1:roi(2)+roi(4)); %XXX correct or exchanged??? a: looks fine
   else
       disp(['could not correct gain and offset maps because camera ' metadata.correctionfile ' was not found']);
   end
end
if isempty(imscmos)
   imphot=(single(imgp)-metadata.offset)*adu2phot;
else
    imphot=imscmos;
end

if correctemmirror && metadata.EMon
    imphot=imphot(:,end:-1:1);
end

end