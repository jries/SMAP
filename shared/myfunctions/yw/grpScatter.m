function ax = grpScatter(varargin)
    % Usage:
    %   grpScatter(ax,x,y,grp,sz,c,varargin)
    %
    if isa(varargin{1}, 'matlab.graphics.axis.Axes')
        ax = varargin{1};
        varargin(1) = [];
    else
        fig = figure;
        ax = axes(ax);
    end
    nargin = length(varargin);
    if nargin<3
        error('grpScatter: data missing or group not specified.');
    else
        x = varargin{1};
        y = varargin{2};
        grp = varargin{3};
        varargin(1:3) = [];
    end
    
    if nargin>3
        sz = varargin{1};
        varargin(1) = [];
    end
    
    if nargin>4
        c = varargin{1};
        varargin(1) = [];
    end
    
    grpID = unique(grp);
    for k = 1:length(grpID)
        lKept = grp==grpID(k);
        switch nargin
            case 3
                scatter(ax, x(lKept), y(lKept))
            case 4
                scatter(ax, x(lKept), y(lKept), sz)
            case 5
                scatter(ax, x(lKept), y(lKept), sz, c(lKept))
            otherwise
                scatter(ax, x(lKept), y(lKept), sz, c(lKept), varargin{:})
        end
        
    hold(ax, 'on')
    end
    hold(ax, 'off')
end