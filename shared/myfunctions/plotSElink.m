function plotSElink(varargin)
% axis (optional), x, y, sitenumber, SE-object, plot-parameters

if isa(varargin{1},'matlab.graphics.axis.Axes')
    startind=2;
    ax=varargin{1};
else
    ax=gca;
    startind=1;
end

xdat0=varargin{startind}(:);
ydat0=varargin{startind+1}(:);
sitenumber0=varargin{startind+2}(:);
[xdat,sortind]=sort(xdat0);
ydat=ydat0(sortind);
sitenumber=sitenumber0(sortind);

SE=varargin{startind+3};
plotpar=varargin{startind+4};
plothandle=plot(ax,xdat,ydat,plotpar);

parent=ax.Parent;
while ~isa(parent,'matlab.ui.Figure')
    parent=parent.Parent;
end
    

dc=datacursormode(parent);
dc.Enable='on';
fndat=dc.UpdateFcn;
% oldupdate=@dc.UpdateFcn;
try
if length(fndat)>1
xdat=[fndat{2}(:),xdat(:)];
ydat=[fndat{3}(:),ydat(:)];
sitenumber=[fndat{4}(:),sitenumber(:)];
end
dc.UpdateFcn={@updatecursor,xdat(:),ydat(:),sitenumber(:),SE};
catch err
    err
end

% bo=brush(ax.Parent);
% bo.ActionPostCallback={@updatebrush,xdat(:),ydat(:),sitenumber(:),SE};
end

% function txt=updatebrush(a,dc,xdat,ydat,sitenumber,SE)
% 
% end
function txt=updatecursor(a,dc,xdat,ydat,sitenumber,SE)
%      if plothandle==dc.Target
    pos=dc.Position;
    ind=find(pos(1)==xdat&pos(2)==ydat);
    if ~isempty(ind)
        site=sitenumber(ind(1));
        txt{1}=['Site ID: ' num2str(site)];
        SE.processors.preview.updatesite(site);
        
    else
        txt={'none'};
%          txto=oldupdate(a,dc)
     end
%     txto=oldupdate(a,dc)

end


