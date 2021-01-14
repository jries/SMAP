function out=myistrue(varargin)
% shorter writing for empty values, fields
% istrue(in): false if isempty(in)
% istrue(p,field): true if field exist and is true
if nargin==1
    if isempty(varargin{1})
        out=false;
    else
        out=logical(varargin{1});
    end
elseif nargin==2
    if isfield(varargin{1},varargin{2}) && ~isempty(varargin{1}.(varargin{2})) && varargin{1}.(varargin{2})
        out=true;
    else
        out=false;
    end
else
    warning('istrue called with wrong number of input parameters')
end
    