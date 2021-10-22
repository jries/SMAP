function plotFreeRot(obj, varargin)
%% PLOTFREEROT Image rotation triggered by the rotatoin of points visulization
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
    
    ax = nexttile(t);           % handle to the scatter plot
    subViz = nexttile(t);       % handle to the image
    disable = uicontrol(ax.Parent.Parent ,'Style','checkbox','String','Disable','Value',1,'Position',[10 360 70 30],'Callback',{@edit_callBack, obj,'disable','value'});
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
        uip = uipanel('Parent',ax.Parent.Parent,'Position',[.01 .01 .98 .2]);
        bg = uibuttongroup('Parent',uip,'Position',[.8 .01 .2 .98]);
        uicontrol(bg ,'Style','text','String','Pixelated','Position',[1 42 60 20]);
        data_sel = uicontrol(bg ,'Style','radiobutton','String','Data','Position',[1 21 60 20],'Callback',{@edit_callBack, obj,'target2Render','string'});
        model_sel = uicontrol(bg ,'Style','radiobutton','String','Model','Position',[1 1 60 20],'Callback',{@edit_callBack, obj,'target2Render','string'});

        section = uicontrol(uip ,'Style','edit','String',50,'Position',[1 26 60 25],'Callback',{@edit_callBack, obj,'section','string'});
        playback = uicontrol(uip ,'Style','pushbutton','String','Playback','Position',[200 26 60 25],'Callback',{@playback_callBack,obj,ax,locs,results,bg, section});
        timeLine = uicontrol(uip ,'Style','slider','Position',[260 26 120 25],'Callback',{@timeLine_callBack,obj,ax,locs, results,bg, section});
        
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
            section.String = obj.display.plotFreeRot.section;
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
        
        locsViz = pointViz(ax, locs, obj);
        obj.setTemp('subViz', subViz);
        obj.setTemp('axViz', ax);
        
        h.ActionPostCallback = {@rotate_callBack, obj, results, bg, section};
        axis(subViz,'equal')
        xlabel(ax, 'Aligned X (nm)')
        ylabel(ax, 'Aligned Y (nm)')
        zlabel(ax, 'Aligned Z (nm)')
        hold(subViz, 'on');
        plot(subViz,  locsViz.xnm+obj.roiSize/2,locsViz.ynm+obj.roiSize/2,' or', 'MarkerEdgeColor','w','MarkerFaceColor','r')
        set(subViz,'YDir','reverse')
        set(ax,'YDir','reverse')
        hold(subViz,'off')
        rotate_callBack([],[],obj,results,bg, section);          % update the image visulization

        % Buttons for the rotation
        uicontrol(uip ,'Style','pushbutton','String','>','Position',[1 1 60 25],'Callback',{@azimuthalRot_callBack, ax, 5, obj,results,bg, section});
        uicontrol(uip ,'Style','pushbutton','String','>>','Position',[30 1 60 25],'Callback',{@azimuthalRot_callBack, ax, 20, obj,results,bg, section});
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
    obj.setTemp('locsViz', locsViz);
    obj.setTemp('modViz', modViz);
end
function rotate_callBack(a,b,obj, results,bg, section)
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
            obj.rotCoordNMkImg(subViz, modViz, locsViz, rotVizAlt, pixelSize, 'Model', str2double(section.String), results.lutLocs);
        case 'Data'
            obj.rotCoordNMkImg(subViz, modViz, locsViz, rotVizAlt, 2, 'Data', str2double(section.String), results.lutLocs);
    end
    set(subViz,'YDir','normal')
    title(subViz,'2D Projection of the left panel')
end

function playback_callBack(a,b,obj,ax,locs, results,bg, section)
    optimHistory = obj.fitInfo.optimHistory;
    numOfIter = size(optimHistory,1);
    for k = 1:numOfIter
        snapshot(obj,ax,locs, k, optimHistory, results,bg, section);
    end
end

function timeLine_callBack(a,b,obj,ax,locs, results,bg, section)
    optimHistory = obj.fitInfo.optimHistory;
    numOfIter = size(optimHistory,1);
    idx = round(numOfIter*a.Value);
    snapshot(obj,ax,locs, idx, optimHistory, results,bg, section)
end

function snapshot(obj,ax,locs, idx, optimHistory, results,bg, section)
lFix = obj.allParsArg.fix;
obj.allParsArg.value(~lFix) = optimHistory(idx,:);
%         w = waitforbuttonpress;
pointViz(ax, locs, obj);
rotate_callBack([],[],obj, results,bg, section);
drawnow
end
%% 
% 
function azimuthalRot_callBack(a,b,ax,dRot,obj,results,bg, section)
    rotVizAlt = get(ax, 'View');
    set(ax, 'View', rotVizAlt+[dRot 0]);
    rotate_callBack([],[],obj,results,bg, section);
end

function edit_callBack(a,b,obj, fn,field)
switch field
    case 'value'
        obj.display.plotFreeRot.(fn) = a.Value;
    otherwise
        obj.display.plotFreeRot.(fn) = a.String;
end
end