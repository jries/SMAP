function luto=lutinvert(lut)
s=size(lut);
luto=lut;
for k=1:s(1)
    luto(k,:)=sum(lut(k,:))-lut(k,:);
end
end