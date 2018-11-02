function indg=brightestinregion(x,y,int,R)
%%

[xs,xind]=sort(x);
ys=y(xind);
ints=int(xind);
[~,iinds]=sort(ints,'descend');

for k=1:length(ints)

    indh=iinds(k);
    if ints(indh)==0
        l=l+1;
        continue
    end
    xh=xs(indh);
    l=1;
    while l<=length(iinds)-indh && xs(indh+l)<xh+R
        if ints(indh+l)==0
            l=l+1;
            continue
        end
        d2=(xs(indh+l)-xh).^2+(ys(indh+l)-ys(indh)).^2;
        if d2<R^2 
            if ints(indh)>=ints(indh+l)
                ints(indh+l)=0;
            else
                ints(indh)=0;
            end
        end
        l=l+1;
    end
    l=1;
    while l<indh&&xs(indh-l)>xh-R
        if ints(indh-l)==0
            l=l+1;
            continue
        end
        d2=(xs(indh-l)-xh).^2+(ys(indh-l)-ys(indh)).^2;
        if d2<R^2 
            if ints(indh)>=ints(indh-l)
                ints(indh-l)=0;
            else
                ints(indh)=0;
            end
        end
        l=l+1;
    end
end
indg=true(size(int));
indg(xind)=ints>0;
