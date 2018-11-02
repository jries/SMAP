function out=mynum2cellstring(in)
out={};
out=in;
for k=1:length(in(:))
    thiselement=in{k};
    if iscell(thiselement)
        thiselement=thiselement{1};
    end
   
    if isnumeric(thiselement)
        out{k}=num2str(thiselement);
    else
        out{k}=thiselement;
    end
end