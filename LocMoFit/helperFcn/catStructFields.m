function S = catStructFields(dim, varargin)
F = cellfun(@fieldnames,varargin,'uni',0);
assert(isequal(F{:}),'All structures must have the same field names.')
T = [varargin{:}];
S = struct();
F = F{1};
for k = 1:numel(F);
    S.(F{k}) = cat(dim,T.(F{k}));
end
end