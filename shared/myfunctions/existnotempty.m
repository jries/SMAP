function out=existnotempty(st,varargin)
out=true;
if isempty(st)
    out=false;
    return
end
for k=1:length(varargin)
    if isempty(st) || ~myisfield(st,varargin{k})
        out=false;
        return
    end
    st=st.(varargin{k});
end
if isempty(st)
    out=false;
end
end