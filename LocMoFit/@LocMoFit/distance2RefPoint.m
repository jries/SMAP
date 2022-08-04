function d = distance2RefPoint(obj, refPoint,varargin)
% :func:`distance2RefPoint` calculates for each localizations the d
%
% Usage:
%   obj.distance2RefPoint(refPoint)
%   obj.distance2RefPoint(refPoint, locs)
% Args:
%   refPoint: the coordinate of the reference point

locs = obj.locs;
if ~isempty(varargin)
    len = length(varargin);
    for k = 1:len
        switch k
            case 1
                locs = varargin{1};
        end
    end
end

end