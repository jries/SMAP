function plotFreeRot(obj, varargin)
% :func:`plotFreeRot` creates the GUI for a user to explore data in 3D with
% free rotations.
% 
% Uasage:
%   obj.plotFreeRot(ax, locs)
%   obj.plotFreeRot(ax, locs)

% PLOTFREEROT Image rotation triggered by the rotatoin of points visulization
% Different layers should have different colours. The lut should be based on 
% what user defined or based on the setting in SMAP.
    % allow user to define the axes
    
    if isa(varargin{1},'matlab.graphics.axis.Axes')
            ax = varargin{1};
            if isa(ax.Parent,'matlab.ui.container.Tab')
                t = tiledlayout(ax.Parent,1,2);
            elseif isa(ax.Parent, 'matlab.graphics.layout.TiledChartLayout')
                t = tiledlayout(ax.Parent.Parent,1,2);
            end
            locs = varargin{2};
            varargin(1:2)=[];
        else
            locs = varargin{1};
            varargin(1)=[];
            fig = figure;
            t = tiledlayout(fig,1,2);
    end
    t.Position = [0.1300 0.2100 0.7750 0.7150];
    ax = nexttile(t);           % handle to the scatter plot
    subViz = nexttile(t);       % handle to the image
    disable = uicontrol(ax.Parent.Parent ,'Style','checkbox','String','Disable','Value',0,'Position',[10 360 70 30],'Callback',{@edit_callBack, obj,'disable','value'});
    if isfield(obj.display, 'plotFreeRot') && isfield(obj.display.plotFreeRot, 'disable')
        disable.Value = obj.display.plotFreeRot.disable;
    end
    
    if disable.Value == 0
        % default settings
        p = inputParser;
        p.addParameter('lutLocs',mymakelut)
        p.addParameter('sigma',12)
        p.addParameter('pixelSize',3)
        p.parse(varargin{:});
        results = p.Results;
        
        % Control panel
        uip = uipanel('Title','Control','Parent',ax.Parent.Parent,'Position',[.01 .01 .98 .2]);
        bg = uibuttongroup('Parent',uip,'Position',[.8 .01 .2 .98]);
        uicontrol(bg ,'Style','text','String','Image soure','Position',[1 42 70 20]);
        guihandles = [];
        data_sel = uicontrol(bg ,'Style','radiobutton','String','Data','Position',[1 21 60 20],'Callback',{@edit_callBack, obj,'target2Render','string'});
        model_sel = uicontrol(bg ,'Style','radiobutton','String','Model','Position',[1 1 60 20],'Callback',{@edit_callBack, obj,'target2Render','string'});

        guihandles.t_section = uicontrol(uip ,...
            'Style','checkbox',...
            'String','Cross-section y +/-',...
            'value',0,...
            'Position',[1 0 1.5 1],...
            'Callback',{@edit_callBack, obj,'t_section','value'});
        guihandles.section = uicontrol(uip ,...
            'Style','edit',...
            'String',50,...
            'Position',[2.2 0 0.4 1],...
            'Callback',{@edit_callBack, obj,'section','string'});
        
        guihandles.t_playback = uicontrol(uip ,'Style','text','String','Optimization playback','Position',[2.6 1 1.2 1]);
        if (~isfield(obj.fitInfo,'optimHistory'))||isempty(obj.fitInfo.optimHistory)
            lEnable = 'off';
        else
            lEnable = 'on';
        end
        guihandles.playback = uicontrol(uip ,'Style','pushbutton','String','|>','Position',[3.8 1 0.3 1],'Callback',{@playback_callBack,obj,ax,locs,results,bg, guihandles.section},'Enable',lEnable);
        guihandles.timeLine = uicontrol(uip ,'Style','slider','Position',[4.1 1 0.9 1],'Callback',{@timeLine_callBack,obj,ax,locs, results,bg, guihandles.section, guihandles.section},'Enable',lEnable);
        
        if isfield(obj.display.plotFreeRot, 'target2Render')
            switch obj.display.plotFreeRot.target2Render
                case 'Model'
                    model_sel.Value = 1;
                    data_sel.Value = 0;
                    bg.SelectedObject = model_sel;
                case 'Data'
                    data_sel.Value = 1;
                    model_sel.Value = 0;
                    bg.SelectedObject = data_sel;
            end
        end

        if isfield(obj.display.plotFreeRot, 'section')
            guihandles.section.String = obj.display.plotFreeRot.section;
        end

        if isfield(obj.display.plotFreeRot, 't_section')
            guihandles.t_section.Value = obj.display.plotFreeRot.t_section;
        end
        
        % Point visualization of the fit
        h = rotate3d(ax);
        h.Enable = 'on';
        obj.setTemp('axVizRot',h)
%         rotVizOri = get(ax, 'View');
%         obj.setTemp('rotVizOri', rotVizOri);
        if isfield(obj.display.plotFreeRot, 'View')
            subViz.View = obj.display.plotFreeRot.View;
            ax.View = obj.display.plotFreeRot.View;
        else
            subViz.View = [0 0];
            ax.View = [0 0];
        end
        
        % Create the scatter plot
        locsViz = pointViz(ax, locs, obj);
        obj.setTemp('subViz', subViz);
        obj.setTemp('axViz', ax);
        
        h.ActionPostCallback = {@rotate_callBack, obj, results, bg, guihandles.section, guihandles.t_section};
        axis(subViz,'equal')
        xlabel(ax, 'Aligned X (nm)')
        ylabel(ax, 'Aligned Y (nm)')
        zlabel(ax, 'Aligned Z (nm)')
        hold(subViz, 'on');
        plot(subViz,  locsViz.xnm+obj.roiSize/2,locsViz.ynm+obj.roiSize/2,' or', 'MarkerEdgeColor','w','MarkerFaceColor','r')
        
        hold(subViz,'off')
        rotate_callBack([],[],obj,results,bg, guihandles.section, guihandles.t_section);          % update the image visulization

        % Buttons for the rotation
        guihandles.t_zrot = uicontrol(uip ,'Style','text','String','Rotation about z','Position',[1 1 1 1]);
        guihandles.zrot_slow = uicontrol(uip ,'Style','pushbutton','String','>','Position',[1.9 1 0.3 0.8],'Callback',{@azimuthalRot_callBack, ax, 5, obj,results,bg, guihandles.section});
        guihandles.zrot_fast = uicontrol(uip ,'Style','pushbutton','String','>>','Position',[2.2 1 0.3 0.8],'Callback',{@azimuthalRot_callBack, ax, 20, obj,results,bg, guihandles.section});
        guiStyle(guihandles, fieldnames(guihandles));
     end
end
%% 
% 
function locsViz = pointViz(ax, locs, obj)
    cla(ax)
    [~,modViz] = obj.plot(ax, locs,'plotType','point','modelSamplingFactor',0.3); % get point type visualization
%     modViz = obj.getLayerPoint(0.75); % get point type visualization
    axes(ax)
    legend({'Model','Data'})
    lPars = obj.exportPars(1,'lPar');
    lPars.variation = 0;
    moveModel = 0;
    if moveModel
        modViz_locFormat.xnm = modViz{1}.x;
        modViz_locFormat.ynm = modViz{1}.y;
        modViz_locFormat.znm = modViz{1}.z;
        modViz_locFormat = obj.locsHandler(modViz_locFormat,lPars,1, 'order_transform','RT');
        modViz{1}.x = modViz_locFormat.xnm;
        modViz{1}.y = modViz_locFormat.ynm;
        modViz{1}.z = modViz_locFormat.znm;
        locsViz = locs;
    else
        locsViz = obj.locsHandler(locs,lPars,1);
    end
    locs.ynm = -locs.ynm;
    obj.setTemp('locsViz', locsViz);
    obj.setTemp('modViz', modViz);
end
function rotate_callBack(a,b,obj, results,bg, section,t_section)
    pixelSize = obj.model{1}.pixelSize;
    pointImg_rotDiff = [0 -90];
    subViz = obj.getTemp('subViz');
    axViz = obj.getTemp('axViz');
    obj.display.plotFreeRot.View = get(axViz, 'View');
    rotVizAlt = get(axViz, 'View')+pointImg_rotDiff;    % this makes sure the point view and the image view will have the same orientation in 3D
    modViz = obj.getTemp('modViz');
    locsViz = obj.getTemp('locsViz');
    
    switch bg.SelectedObject.String
        % either show the data or the model as an image
        case 'Model'
            obj.rotCoordNMkImg(subViz, modViz, locsViz, rotVizAlt, pixelSize, 'Model', str2double(section.String)./t_section.Value, results.lutLocs);
        case 'Data'
            obj.rotCoordNMkImg(subViz, modViz, locsViz, rotVizAlt, 2, 'Data', str2double(section.String)./t_section.Value, results.lutLocs);
    end
    title(subViz,'2D Projection of the left panel')
    set(subViz,'YDir','normal')
end

function playback_callBack(a,b,obj,ax,locs, results,bg, section, t_section)
    optimHistory = obj.fitInfo.optimHistory;
    numOfIter = size(optimHistory,1);
    for k = 1:numOfIter
        snapshot(obj,ax,locs, k, optimHistory, results,bg, section, t_section);
    end
end

function timeLine_callBack(a,b,obj,ax,locs, results,bg, section, t_section)
    optimHistory = obj.fitInfo.optimHistory;
    numOfIter = size(optimHistory,1);
    idx = round(numOfIter*a.Value);
    snapshot(obj,ax,locs, idx, optimHistory, results,bg, section, t_section)
end

function snapshot(obj,ax,locs, idx, optimHistory, results,bg, section, t_section)
lFix = obj.allParsArg.fix;
obj.allParsArg.value(~lFix) = optimHistory(idx,:);
%         w = waitforbuttonpress;
pointViz(ax, locs, obj);
rotate_callBack([],[],obj, results,bg, section, t_section);
drawnow
end
%% 
% 
function azimuthalRot_callBack(a,b,ax,dRot,obj,results,bg, section, t_section)
    rotVizAlt = get(ax, 'View');
    set(ax, 'View', rotVizAlt+[dRot 0]);
    rotate_callBack([],[],obj,results,bg, section, t_section);
end

function edit_callBack(a,b,obj, fn,field)
switch field
    case 'value'
        obj.display.plotFreeRot.(fn) = a.Value;
    otherwise
        obj.display.plotFreeRot.(fn) = a.String;
end
end