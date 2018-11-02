function out=myisfield(object,name)
if isstruct(object)
    out=isfield(object,name);
else
    out=isprop(object,name);
end
