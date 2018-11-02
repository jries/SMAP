function out=myvalidatestring(x,fields)
if ischar(x)
    out=validatestring(x,fields);
elseif iscell(x)
    for k=1:length(x)
        out=validatestring(x{k},fields);
    end
elseif isnumeric(x)  ||islogical(x)
    
    out=1;
end
end

%handle also cell x
