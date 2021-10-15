function fig = plotFixRot(obj, varargin)
%% PLOTFIXROT Plot 2D projections of the 3D fit at fixed orientations
% Different layers should have different colours. The lut should be based on
% what user defined or based on the setting in SMAP.
% allow user to define the axes
if isa(varargin{1},'matlab.ui.container.Tab')
    fig = varargin{1};
    locs = varargin{2};
    varargin(1:2)=[];
    ax = axes(fig);
    ax.Visible = 'off';
else
    locs = varargin{1};
    varargin(1)=[];
    fig = figure;
end

disable = uicontrol(fig ,'Style','checkbox','String','Disable','Value',1,'Position',[10 360 70 30],'Callback',{@edit_callBack, obj,'disable','value'});
if isfield(obj.display, 'plotFixRot') && isfield(obj.display.plotFixRot, 'disable')
    disable.Value = obj.display.plotFixRot.disable;
end

if disable.Value == 0
    % default settings
    p = inputParser;
    p.addParameter('lutLocs',mymakelut)
    p.addParameter('sigma',12)
    p.addParameter('pixelSize',3)
    p.parse(varargin{:});
    results = p.Results;
    
    %% define layout
    % for 4 different orientations
    t = tiledlayout(fig, 1,4);
    subViz1 = nexttile(t);
    subViz2 = nexttile(t);
    subViz3 = nexttile(t);
    subViz4 = nexttile(t);
    
    % Point visualization of the fit
    for j = obj.numOfModel:-1:1
        variation_ori(j) = obj.getVariable(['m' num2str(j) '.lPar.variation']);
        obj.setParArg(['m' num2str(j) '.lPar.variation'], 'value',0);
    end
    modViz = obj.getModPoint(0.3); % get point type visualization
    locsViz = obj.locsHandler(locs,obj.exportPars(1,'lPar'),1);
    obj.rotCoordNMkImg(subViz1, modViz, locsViz, [0 -90], 2, 'Data', 30, results.lutLocs)
    obj.rotCoordNMkImg(subViz2, modViz, locsViz, [45 -90], 2, 'Data', 30, results.lutLocs)
    obj.rotCoordNMkImg(subViz3, modViz, locsViz, [90 -90], 2, 'Data', 30, results.lutLocs)
    obj.rotCoordNMkImg(subViz4, modViz, locsViz, [135 -90], 2, 'Data', 30, results.lutLocs)
    set(subViz1,'YDir','normal')
    axis(subViz1,'equal')
    set(subViz2,'YDir','normal')
    axis(subViz2,'equal')
    set(subViz3,'YDir','normal')
    axis(subViz3,'equal')
    set(subViz4,'YDir','normal')
    axis(subViz4,'equal')
    for j = obj.numOfModel:-1:1
        obj.setParArg(['m' num2str(j) '.lPar.variation'], 'value',variation_ori(j));
    end
end
end

function edit_callBack(a,b,obj, fn,field)
switch field
    case 'value'
        obj.display.plotFixRot.(fn) = a.Value;
    otherwise
        obj.display.plotFixRot.(fn) = a.String;
end
end