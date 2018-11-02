%
function outim=cutoutim(in,posin)
%position=[left right bottom top]

s=size(in);
if length(s)>2
    for k=1:length(s)
        im3{k}=cutoutim2D(in(:,:,k),posin);
    end
    s2=size(im3{1});
    outim=zeros(s2(1),s2(2),s(3));
    for k=1:s(3)
        outim(:,:,k)=im3{k};
    end
else
    outim=cutoutim2D(in,posin);
end


function outim=cutoutim2D(in,posin)
position=round(posin);
extralines=[0 0 0 0];
sin=size(in);
if position(1)<1
    extralines(1)=-position(1)+1;
    position(1)=1;
end
if position(3)<1
    extralines(3)=-position(3)+1;
    position(3)=1;
end
if position(2)>sin(1)
    extralines(2)=position(2)-sin(1);
    position(2)=sin(1);
end
if position(4)>sin(2)
    extralines(4)=position(4)-sin(2);
    position(4)=sin(2);
end

outim=in(position(1):position(2),position(3):position(4));

so=size(outim);
outim=[zeros(extralines(1),so(2)); outim];
so=size(outim);
outim=[outim;zeros(extralines(2),so(2))];
so=size(outim);
outim=[zeros(so(1),extralines(3)), outim];
so=size(outim);
outim=[outim,zeros(so(1),extralines(4))];
