function badout=tooclose(x,y,mind)

mask=false(size(x));
for k=1:length(x)
    if mask(k)
        continue
    end
    bad=(x-x(k)).^2+(y-y(k)).^2<mind^2;
    if sum(bad)>1
        mask(bad)=true;
    end
end

badout=mask;
end