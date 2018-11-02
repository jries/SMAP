function to=mytextsplit(txt,width,tol)
if nargin<3
    tol=8;
end
if nargin<2
    width=100;
end

ind=1;
indl=1;
to='';
while 1
    searchrange=max(1,ind+width-tol):min(ind+width+tol,length(txt));
    if isempty(searchrange)
        to=[to 10 txt(ind:end)];
        break;
    end
    split=find(txt(searchrange)<=32,1,'first')+searchrange(1)-1;
    if isempty(split)
        split=ind+width;
        if split>length(txt)
            to=[to 10 txt(ind:end)];
            break
        end
    end
    txt(ind:split);
    to=[to 10 txt(ind:split)];
    ind=split+1;
    indl=indl+1;
end
