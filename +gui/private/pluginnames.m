
function out=pluginnames(varargin)
plugs=plugin;
for k=1:length(varargin)
    plugs=plugs.(varargin{k});
end
if isstruct(plugs)
out=fieldnames(plugs);
else
    out=plugs;
end
end