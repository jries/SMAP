function [hiso, hcap, oneFrame] = patchPlot(varargin)
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
p.addParameter('azi_light',[]);
p.addParameter('ele_light',[]);
p.addParameter('azi_view',[]);
p.addParameter('ele_view',[]);
p.addParameter('colormap','winter');
p.addParameter('FaceColor','#C0C0C0');
p.addParameter('isoCutoff',1.1);% for isosurface, at max(V(:))/cutoff
p.parse(varargin{:});
p = p.Results; 
p.sgauss=30;%smoothing of volume
% p.cutoff=1.2; % for isosurface, at max(V(:))/cutoff
p.cutoffdc=[3 3]; % for isosurface, at max(V(:))/cutoff
p.pxSize = 1; % pixel size

p.smooth = false;
if p.smooth
    vs=smooth3(v,'gaussian',[1 1 1]*(round(2*p.sgauss/2)*2+1),p.sgauss);
else
    vs = v;
end
co=max(vs(:))/p.isoCutoff;

fc = [1,.75,.65];
%             fc = 'red';
% [f,v] = isosurface(vs,co);
% [f2,v2] = isosurface(p.v2,co);
[f1,v1] = isosurface(vs,co);
hold(ax, 'on')
hiso = patch(ax, 'Faces',f1,'Vertices',v1 ,'EdgeColor','none','FaceColor',p.FaceColor);
hold(ax, 'off')
isonormals(vs,hiso);
hiso.DiffuseStrength = 0.9;
hiso.SpecularStrength = 0.1;
hiso.SpecularExponent = 15;
%% Enable these two lines to show the 2nd half of the sphere 
axis equal
axis vis3d
% xlim([0 range(my)/p.pxSize])
% ylim([0 range(mx)/p.pxSize])
% zlim([0 range(mz)/p.pxSize])
view(p.azi_view,p.ele_view)
colormap winter
%             xlim([0,size(vs,1)])
h = light(ax);
lightangle(h,p.azi_light,p.ele_light);
lighting gouraud
colormap(p.colormap);
drawnow
oneFrame = getframe(fig);
oneFrame = oneFrame.cdata;
end