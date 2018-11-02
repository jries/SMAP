function out=mywaveletfilteratrous(in,filter)
% inspired by [1]	I. Izeddin, J. Boulanger, V. Racine, C. G. Specht, A. Kechkar, D. Nair, A. Triller, D. Choquet, M. Dahan, and J. B. Sibarita, ?Wavelet analysis for single molecule localization microscopy,? Opt Express, vol. 20, no. 3, pp. 2081?2095, Jan. 2012.
persistent H0 H1 H2 g1 g2
if isempty(H0)
    H0=3/8;
    H1=1/4;
    H2=1/16;
    g1=[H2,H1,H0,H1,H2];
    g2=[H2,0,H1,0,H0,0,H1,0,H2];
end

if 0
    co=2*mean(in(:));
   
    in(in>co)=co;  
end

V1=conv2(conv2(in,g1','same'),g1,'same');
V2=conv2(conv2(V1,g2','same'),g2,'same');
if filter
    out=V1-V2;
else
    out=V2;
end


% figure(88);
% h=fspecial('gaussian',5,1);
% imagesc(horzcat(in,V1,V2,V1-V2,out,filter2(h,out)))
end