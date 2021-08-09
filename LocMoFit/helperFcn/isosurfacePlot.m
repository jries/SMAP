function [hiso, hcap, oneFrame] = isosurfacePlot(varargin)
if isa(varargin{1},'matlab.graphics.axis.Axes')
    ax = varargin{1};
    fig = ax.Parent;
    varargin(1) = [];
else
    fig = figure;
    ax = axes(fig);
end
if length(size(varargin{1}))==3&&ismatrix(varargin{1}(:,:,1))
    v = varargin{1};
    varargin(1) = [];
else
    warning('Input should be a 3D matrix.')
    return
end
p = inputParser;
p.addParameter('v2',[]);
p.parse(varargin{:});
p = p.Results; 
p.sgauss=20;%smoothing of volume
p.cutoff=4; % for isosurface, at max(V(:))/cutoff
p.cutoffdc=[3 3]; % for isosurface, at max(V(:))/cutoff
p.pxSize = 5; % pixel size
p.smooth = false;
if p.smooth
    vs=smooth3(v,'gaussian',[1 1 1]*(round(2*p.sgauss/2)*2+1),p.sgauss);
else
    vs = v;
end
co=max(vs(:))/p.cutoff;

fc = [1,.75,.65];
%             fc = 'red';
% [f,v] = isosurface(vs,co);
% [f2,v2] = isosurface(p.v2,co);
[f1,v1] = isosurface(vs(end/2:end,:,:),co);
[f2,v2] = isosurface(p.v2(1:end/2,:,:),co);
v1(:,2) = v1(:,2)+size(vs,2)/2-2;
v2(:,2) = v2(:,2);
hiso = patch(ax, 'Faces',f1,'Vertices',v1 ,'EdgeColor','none','FaceColor','#C0C0C0');
isonormals(vs,hiso);

%% Enable these two lines to show the 2nd half of the sphere 
hiso2 = patch(ax, 'Faces',f2,'Vertices',v2,'EdgeColor','none','FaceColor','#C0C0C0','FaceAlpha',.3);
isonormals(p.v2,hiso2);

%%
hold on
[f,v,c]= isocaps(vs(end/2:end,:,:),co);
v(:,2) = v(:,2)+size(vs,2)/2-2;
hcap = patch(ax, 'Faces',f,'Vertices',v,'FaceVertexCData',c,'FaceColor','interp','EdgeColor','none');
% ,'FaceColor','interp'
axis equal
% xlim([0 range(my)/p.pxSize])
% ylim([0 range(mx)/p.pxSize])
% zlim([0 range(mz)/p.pxSize])
view(-35,30)
colormap winter
%             xlim([0,size(vs,1)])
lightangle(45,30);
lighting gouraud
colormap winter;
drawnow
oneFrame = getframe(fig);
oneFrame = oneFrame.cdata;
end