function z=zfromSx2_Ss2(sx2sy,sx,sy)
sx2sy2=sx.^2-sy.^2;
z=feval(sx2sy,sx2sy2);
z(sy==0)=-inf;
                
                
                
%                 
%  sxr=round((sx-zmat.smin)/zmat.ds)+1;
% syr=round((sy-zmat.smin)/zmat.ds)+1;
% 
% ns=length(zmat.s);
% goodind=~(sxr<1|sxr>ns|syr<1|syr>ns);
% 
% sind=sub2ind(size(zmat.z),sxr(goodind),syr(goodind));
% z=zeros(size(sx),'single');
% z(goodind)=zmat.z(sind);
% z(~goodind)=-2;
% d=zeros(size(sx),'single');
% d(goodind)=zmat.d(sind);
% d(~goodind)=inf;
