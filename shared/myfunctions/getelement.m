function out=getelement(in,varargin)
if iscell(in)
    out=in{varargin{:}};
else
    out=in(varargin{:});
end