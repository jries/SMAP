function [rs,norm]=radialsum(img)
%from FRC, Rieger
    %FRCresolution calculates FRC resolution
    %     FRC resolution implementation based on the matlab code provided with: 
    %[1]	R. P. J. Nieuwenhuizen, K. A. Lidke, M. Bates, D. L. Puig, D. Gr?nwald,
    %S. Stallinga, and B. Rieger, ?Measuring image resolution in optical nanoscopy.,? 
    %Nat Methods, vol. 10, no. 6, pp. 557?562, Jun. 2013.
s=size(img);
center=floor((s+1)/2);
rs=zeros(ceil(s(1)/2)+1,1);
norm=zeros(ceil(s(1)/2)+1,1);
for k=1:s(1)
    for l=1:s(2)
        d=sqrt((k-center(1)).^2+(l-center(2)).^2);
        ind=round(d)+1;
        if ind<=length(rs)
        rs(ind)=rs(ind)+img(k,l);
        norm(ind)=norm(ind)+1;
        end
    end
end
end