%NMS 2D block after  Neubeck and Gool, Algorithm 4

%reprogram in c: much faster hopefully. Now 7ms for 100x100
function maximaout=NMS2DBlockCcall(imin,n)
%maxima=[x,y,intensity]

simg=size(imin);
maxmax=floor(simg(1)*simg(2)/(2*n)^2)*4;
% asdf
       
        
        
maxima= NMS2DBlockC(single(imin),int16(n),uint32(maxmax));
% size(maxima)
% sum(maxima)
f=find(maxima(:,1)==0,1,'first');
maximaout=maxima(1:f-1,:);


% maximaout=maxima(1:maxind-1,:);