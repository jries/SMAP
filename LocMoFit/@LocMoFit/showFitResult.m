function showFitResult(obj, varargin)
    if isa(varargin{1},'matlab.graphics.axis.Axes')
            vis = varargin{1};
            locs = varargin{2};
            varargin(1:2)=[];
        else
            locs = varargin{1};
            varargin(1)=[];
            fig = figure;
            vis = axes(fig);
    end
    newLocs = obj.locsHandler(locs, obj.exportPars(1,'lPar'), [], 'onlyLocpre', true);
    
    v_pixelSize= obj.getTemp('pixelSize_showFitResult');
    if isempty(v_pixelSize)
        v_pixelSize = obj.model{1}.pixelSize;
        obj.setTemp('pixelSize_showFitResult',v_pixelSize);
    end
    v_moveModel = obj.getTemp('moveModel_showFitResult');
    if isempty(v_moveModel)
        v_moveModel = false;
        obj.setTemp('moveModel_showFitResult',v_moveModel);
    end
    
    ax = obj.plot(newLocs,'plotType','image','pixelSize',obj.getTemp('pixelSize_showFitResult'),'movModel',obj.getTemp('moveModel_showFitResult'));
    ax.Position = [0.1300 0.2100 0.7750 0.7150];
    
    h = ax.Parent;
    ax.Parent = vis.Parent;
    delete(vis);
    close(h);
    
    % Control panel
%     uip = uipanel('Title','Control','FontSize',12);
    pard.uip = uipanel('Parent',ax.Parent,'Title','Control','Position',[.01 .01 .98 .2]);
    
    
    if (~isfield(obj.fitInfo,'optimHistory'))||isempty(obj.fitInfo.optimHistory)
        lEnable = 'off';
    else
        lEnable = 'on';
    end
    
    % the unit here is 'line'
    pard.t_playback = uicontrol(pard.uip,...
        'Style','text',...
        'String','Optimization playback',...
        'Position',[1 1 1.5 1]);
    
    pard.playback = uicontrol(pard.uip,...
        'Style','pushbutton',...
        'String','|>',...
        'Position',[2.2 1 0.2 1],...
        'Callback',{@playback_callBack,obj,ax,locs},...
        'Enable', lEnable);
    
    pard.timeLine = uicontrol(pard.uip,...
        'Style','slider',...
        'Position',[2.4 1 1 1],...
        'Callback',{@timeLine_callBack,obj,ax,locs},...
        'Enable', lEnable);
    
    pard.t_pixelSize = uicontrol(pard.uip ,...
        'Style','text',...
        'String','Pixel size',...
        'Position',[3.7 1 0.5 1],...
        'HorizontalAlignment','left');
    
    pard.pixelSize = uicontrol(pard.uip ,...
        'Style','edit',...
        'Position',[4.2 1 0.5 1],...
        'String',v_pixelSize,...
        'Callback',{@setVizPar,obj,'pixelSize_showFitResult'});
    
    pard.moveModel = uicontrol(pard.uip ,...
        'Style','checkbox',...
        'String','Move model',...
        'Value',v_moveModel,...
        'Position',[5 1 1 1],...
        'Callback',{@setVizPar,obj,'moveModel_showFitResult'});
    
    obj.saveHandles(struct('showFitResult',ax.Parent));
    obj.saveHandles(pard, 'showFitResult');
	
    % convert to absolute positions
    guiStyle(obj.handles, strcat(fieldnames(pard),'_showFitResult'), 'exclude', 'uip_showFitResult');
end

function playback_callBack(a,b,obj,ax,locs)
    optimHistory = obj.fitInfo.optimHistory;
    numOfIter = size(optimHistory,1);
    for k = 1:numOfIter
        snapshot(obj,ax,locs, k, optimHistory);
    end
end

function timeLine_callBack(a,b,obj,ax,locs)
    optimHistory = obj.fitInfo.optimHistory;
    numOfIter = size(optimHistory,1);
    idx = max(round(numOfIter*a.Value),1);
    snapshot(obj,ax,locs, idx, optimHistory)
end

function snapshot(obj,ax,locs, idx, optimHistory)
lFix = obj.allParsArg.fix;
obj.allParsArg.value(~lFix) = optimHistory(idx,:);
%         w = waitforbuttonpress;
newLocs = obj.locsHandler(locs, obj.exportPars(1,'lPar'), [], 'onlyLocpre', true);
obj.plot(ax,newLocs,'plotType','image','pixelSize',obj.getTemp('pixelSize_showFitResult'),'movModel',obj.getTemp('moveModel_showFitResult'));
drawnow
end

function setVizPar(a,b,obj,fn)
    switch a.Style
        case 'edit'
            obj.setTemp(fn, str2double(a.String));
        case 'dropdown'
            obj.setTemp(fn, str2double(a.String(a.Value)));
        case 'checkbox'
            obj.setTemp(fn, a.Value);
    end
end