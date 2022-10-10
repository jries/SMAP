function varargout = spacewin_gui(data,params)
% SPACEWIN_GUI a helper gui for making spatial window objects for
% point data
%
% SPACEWIN_GUI(DATA)    Start the gui using the data contained in DATA.
%                       DATA should be a struct or struct array with fields
%                       x, y, and t corresponding to the x, y coordinates
%                       and time of observation of each point of the
%                       dataset. Each struct in an array of structs is
%                       interpreted as belonging to a different channel of
%                       data. The gui provides options to save the resulting
%                       spatial window struct to the base workspace or to a
%                       MAT-file.
%
% SW = SPACEWIN_GUI(_)  As above, but the spatial window struct is returned
%                       as the output arg SW.
%
% H = SPACEWIN_GUI(_,'FigHandle',true) 
% [H,SW] = SPACEWIN_GUI(_,'FigHandle',true) In these forms, a figure handle
%                                       for the gui is also returned as H.
%
% SPACEWIN_GUI(DATA,'Name',Value) Currently two name-value pairs are
%                                 supported:
%     'PixelSize' (scalar double)   pixel size to use to display the
%                                   reconstructed image, in the same units
%                                   as the x and y data.
%     'Ref' (imref2d)               imref2d object to use for
%                                   reconstructed image.

% Copyright (C) 2022 Thomas R Shaw, Sarah L Veatch
% This file is part of SMLM SPACETIME RESOLUTION
% SMLM SPACETIME RESOLUTION is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% SMLM SPACETIME RESOLUTION is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% You should have received a copy of the GNU General Public License
% along with SMLM SPACETIME RESOLUTION.  If not, see <https://www.gnu.org/licenses/>
arguments
    data (1,:) struct {mustHaveXYFields}
    params.PixelSize (1,1) double = 12
    params.Ref = []
    params.FigHandle (1,1) logical = false;
end

% whether we are returning the window as output arg or giving buttons for saving
show_save_btns = (nargout == 0) || (params.FigHandle);

% set up the figure
fig = figure('Name','Spatial Window gui');
set(fig, 'DefaultUicontrolUnits', 'centimeters');
set(fig, 'DefaultUicontrolHorizontalAlignment', 'center');
fig.SizeChangedFcn = @resize_callback;
fig.Units = 'centimeters';
handles = guidata(fig);
% left panel for buttons etc, right panel just for the axes
handles.rpanel = uipanel(fig, 'Units', 'centimeters');
handles.lpanel = uipanel(fig, 'Units', 'centimeters');

% Set up the axes options how we want them
ax = axes(handles.rpanel, 'Units', 'normalized', 'OuterPosition', [0.01 0.01 .98 .98]);
handles.ax = ax;



% Show the data
handles.data = data;
if isempty(params.Ref)
    handles.ref = default_iref([data(1).x(:),data(1).y(:)],params.PixelSize);
else
    handles.ref = params.Ref;
    params.PixelSize = params.Ref.PixelExtentInWorldX;
end
handles.psize = handles.ref.PixelExtentInWorldX;
handles.I = gausblur(reconstruct(data(1),handles.ref), 15/handles.psize);

handles.im_obj = imshow(handles.I,handles.ref, 'Parent', handles.ax);

% make the axes nicer?
caxis(ax,[0,1]);
hold(ax, 'on');
axis(ax,'tight');
ax.Toolbar.Visible = 'off';
ax.Interactions = [zoomInteraction, regionZoomInteraction];
handles.p = polyshape();
handles.p_obj = [];

% set up the various buttons
handles.addbtn = uicontrol(handles.lpanel, 'Style', 'pushbutton', ...
     'String', 'add polygon', 'Callback', @addmask_callback);
handles.rmbtn = uicontrol(handles.lpanel, 'Style', 'pushbutton', ...
    'String', 'remove polygon', 'Callback', @rmmask_callback);

handles.type_label = uicontrol(handles.lpanel, 'Style', 'text', ...
    'String', 'spatial window type');
handles.type_menu = uicontrol(handles.lpanel, 'Style', 'popupmenu',...
    'String',{'polyshape', 'image', 'polygon'}, 'Value', 2);

handles.cmax_label = uicontrol(handles.lpanel, 'Style', 'text', ...
    'String', 'saturation intensity (cmax)');
handles.cmax_edit = uicontrol(handles.lpanel, 'Style', 'edit', ...
        'String', '1', 'Callback', @updatecmax);
handles.channel_label = uicontrol(handles.lpanel, 'Style', 'text', ...
    'String', 'channel number');

Nchannels = numel(data);
channelnames = cell(1,Nchannels);
for i=1:Nchannels, channelnames{i} = num2str(i); end
handles.channel_menu = uicontrol(handles.lpanel, 'Style', 'popupmenu',...
    'String',channelnames, 'Value', 1, 'Callback', @updateimage);
handles.timerange_label = uicontrol(handles.lpanel, 'Style', 'text', ...
    'String', 'time range (min)');
handles.timerange_min = uicontrol(handles.lpanel, 'Style', 'edit',...
    'String',num2str(min(data(1).t)/60, 2),'Callback', @updateimage);
handles.timerange_max = uicontrol(handles.lpanel, 'Style', 'edit',...
    'String',num2str(max(data(1).t)/60, 2),'Callback', @updateimage);

handles.lpanel_uic_list = {'addbtn', 'rmbtn', 'break', 'type_label',...
    'type_menu', 'cmax_label', 'cmax_edit', 'channel_label',...
    'channel_menu', 'timerange_label', 'timerange_min', 'timerange_max'};

if show_save_btns
    handles.varname_label = uicontrol(handles.lpanel, 'Style', 'text', ...
        'String', 'var name for saving');
    handles.varname_edit = uicontrol(handles.lpanel, 'Style', 'edit', ...
        'String', 'spacewin');
    handles.saveasvarbtn = uicontrol(handles.lpanel, 'Style', 'pushbutton',...
        'String', 'Save to base workspace', 'Callback', @savetoworkspace);

    handles.fname_label = uicontrol(handles.lpanel, 'Style', 'text', ...
        'String', 'file name for saving');
    handles.fname_edit = uicontrol(handles.lpanel, 'Style', 'edit', ...
        'String', 'window.mat');
    handles.saveasfilebtn = uicontrol(handles.lpanel, 'Style', 'pushbutton',...
        'String', 'Save to file', 'Callback', @savetomatfile);

    handles.lpanel_uic_list = [handles.lpanel_uic_list,...
        {'varname_label', 'varname_edit', 'saveasvarbtn',...
        'fname_label', 'fname_edit', 'saveasfilebtn'}];
else
    handles.close_btn = uicontrol(handles.lpanel, 'Style', 'pushbutton',...
        'String', 'Close and return spacewin', 'Callback', @(h,e) set(h, 'UserData', true));
    handles.lpanel_uic_list = [handles.lpanel_uic_list,{'break','close_btn'}];
    % make it so closing the window also returns the spacewin
    fig.CloseRequestFcn = @(h,e) set(getfield(guidata(h), 'close_btn'), 'UserData', true);
end

guidata(fig, handles);

pos = fig.Position;
pos = pos + [0 -1 1 1]*5;
set(fig,'Position',pos); % this triggers resizecallback which sets up the correct uicontrol Positions

% return figure handle if requested, or spatial window otherwise
if nargout > 0
    if params.FigHandle
        varargout = {fig};
    else
        fprintf('Note: command line will be unavailable until you press "close and return"\n');
        waitfor(handles.close_btn, 'UserData', true);
        types = handles.type_menu.String;
        type = types{handles.type_menu.Value};
        spacewin = polyshape_to_spacewin(handles.p, type, handles.ref);
        varargout = {spacewin};
        delete(fig);
    end
end

% HELPER FUNCTIONS (mostly callbacks)
    % for recalculating gui layout after the window is resized
    function resize_callback(obj, ~)
        handles = guidata(obj);
        fullpos = obj.Position;
        % left panel is first 5 cm, full height
        lppos = [0,0,5,fullpos(4)];
        handles.lpanel.Position = lppos;
        % right panel is rest
        rppos = [5,0,fullpos(3)-5,fullpos(4)];
        handles.rpanel.Position = rppos;
        
        % layout the uicontrols, vertically.
        left = .2;
        w = 4.6;
        h = .8;
        dy = 1;
        dy_break = .6;
        top = handles.lpanel.Position(4) - dy;
        y = top;
        % just lay them out vertically, separated by dy, taking up (almost)
        % full width of panel
        ll = handles.lpanel_uic_list;
        for j=1:numel(ll)
            uiname = ll{j};
            switch uiname
                case 'break'
                    y = y - dy_break;
                    continue;
                otherwise
                    uielement = handles.(uiname); % pull out handle corresponding to this name
                    uielement.Position = [left, y, w, h];
                    y = y - dy;
                    if strcmp(uielement.Style, 'text')
                        y = y + dy/2;
                    end
            end
        end
    end

    % handle "add mask" button
    function addmask_callback(obj, ~)
        handles = guidata(obj);
        ax = handles.ax;
        
        P = drawpolygon(ax);
        
        handles.p = union(handles.p, polyshape(P.Position));
        delete(P)
        delete(handles.p_obj)
        handles.p_obj = plot(ax,handles.p, 'FaceColor', lines(1));
        guidata(obj, handles);
    end

    % handle "remove mask" button
    function rmmask_callback(obj, ~)
        handles = guidata(obj);
        ax = handles.ax;
        
        P = drawpolygon(ax);
        
        handles.p = subtract(handles.p, polyshape(P.Position));
        delete(P)
        delete(handles.p_obj)
        handles.p_obj = plot(ax,handles.p, 'FaceColor', lines(1));
        guidata(obj, handles);
    end

    function savetoworkspace(obj, ~)
        handles = guidata(obj);
        
        types = handles.type_menu.String;
        type = types{handles.type_menu.Value};
        spacewin = polyshape_to_spacewin(handles.p, type, handles.ref);
        varname = handles.varname_edit.String;
        
        try
            q = evalin('base', varname); %#ok<NASGU>
            %if you get to here, there is a var varname in base
            answer = questdlg(...
                sprintf('Are you sure you want to overwrite variable %s in base workspace?', varname),...
                'Overwrite?', 'Yes', 'Cancel', 'Cancel');
            if strcmp(answer, 'Yes')
                assignin('base', varname, spacewin);
            end
        catch % no var varname in base, proceed
            assignin('base', varname, spacewin);
        end
    end

    function savetomatfile(obj, ~)
        handles = guidata(obj);
        
        types = handles.type_menu.String;
        type = types{handles.type_menu.Value};
        spacewin = polyshape_to_spacewin(handles.p, type, handles.ref);
        varname = handles.varname_edit.String;
        fname = handles.fname_edit.String;
        s.(varname) = spacewin;
        
        if exist(fname, 'file')
            answer = questdlg(...
                sprintf('Are you sure you want to overwrite file %s?', fname),...
                'Overwrite?', 'Yes', 'Cancel', 'Cancel');
            if strcmp(answer, 'Yes')
                save(fname, '-struct', 's', '-append');
            end
        else % no var varname in base, proceed
            save(fname, '-struct', 's');
        end
    end
end

function mustHaveXYFields(v)
if ~isfield(v, 'x') || ~isfield(v,'y')
    error('Struct input must have fields x and y')
end
end

function mustHaveEqualSize(v,u)
if ~isequal(size(v),size(u))
    error('Inputs must be equal sizes')
end
end

function data = resolvefnamevsxy(fname,x,y)
hasfname = ~isempty(fname);
hasx = ~isempty(x);
if ~xor(hasfname, hasx)
    error('spacewin_gui: expecting exactly one nonempty dataset, specified as filename, x and y arrays, or structure with x,y fields')
end
if hasfname
    s = load(fname);
    data = s.data;
end
if hasx
    data.x = x; data.y = y;
end
end

function updatecmax(obj, ~)
        handles = guidata(obj);
        cmax = str2double(handles.cmax_edit.String);
        if isfinite(cmax)
            caxis(handles.ax,[0,cmax]);
        else
            handles.cmax_edit.String = "1";
        end
end

function updateimage(obj, ~)
    handles = guidata(obj);
    channel = handles.channel_menu.Value;
    tmin = 60*str2double(handles.timerange_min.String);
    tmax = 60*str2double(handles.timerange_max.String);
    keep = handles.data(channel).t >= tmin & handles.data(channel).t <= tmax;
    d.x = handles.data(channel).x(keep);
    d.y = handles.data(channel).y(keep);
    handles.I = gausblur(reconstruct(d,handles.ref), 15/handles.psize);
    %handles.im_obj = imshow(handles.I,handles.ref, 'Parent', handles.ax);
    set(handles.im_obj, 'CData', handles.I);
    drawnow
end
