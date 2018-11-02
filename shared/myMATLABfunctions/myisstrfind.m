function out=myisstrfind(a,b)
    if ~iscell(b)
        out=myisstrfindi(a,b);
    else
        out=false;
        for k=1:length(b)
            if myisstrfind(a,b{k});
                out=true;
                return
            end
        end
    end
end

function out= myisstrfindi(a,b)
out=false;
if iscell(a)
    for k=1:length(a)
        if ~isempty(strfind(a{k},b))
            out=true;
            return
        end
    end
else
    if ~isempty(strfind(a,b))
        out=true;
    end   
end
end