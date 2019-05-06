classdef MarkdownPanel < hgsetget & dynamicprops
    %# MarkdownPanel
    % Control which displays markdown as HTML within a MATLAB control
    %
    % ------
    % This control utilizes the [Showdown javascript library][1] to convert
    % markdown into HTML and then uses MATLAB's own `HTMLBrowserPanel` to
    % display this HTML.
    %
    % It behaves like any other graphics object within MATLAB in that all
    % properties can either be set upon object construction
    %
    %     h = MarkdownPanel('Parent', figure, 'Content', '# Hello World!');
    %
    % Or after object creation using the returned handle
    %
    %     h = MarkdownPanel();
    %     h.Parent = gcf;
    %     set(h, 'Position', [0, 0, 0.5, 0.5])
    %
    % To set the actual Markdown content, use the `Content` property. You
    % can provide *either* a string, or a cell array of strings which will
    % automatically create a multi-line entry
    %
    %     set(h, 'Content', '#Hello World')
    %     set(h, 'Content', {'#Hello World', 'This is a test...'})
    %
    % You can use the `Options` property to modify options that are
    % specific to how Showdown renders the markdown. By default, we use all
    % of the default settings except that we enable support for tables.
    %
    %     h.Options.tables = true;
    %
    % The `Options` property is simply a struct where the fieldnames are
    % the option names and the value is the option value. You can modify
    % this struct to adjust an option.
    %
    %     % Enable support for tasklists
    %     h.Options.taskslists = true
    %
    % A complete list of options can be found in the [Showdown
    % documentation][1]
    %
    % ------
    % **Usage**
    %
    %     panel = MarkdownPanel();
    %
    % **Outputs**
    %
    %   `panel`,    Graphics Object, The graphics handle that can be used
    %               to manipulate the appearance of the control
    %
    % ------
    % **Demo**
    %
    % A demo application has been bundled with this code to show how to use
    % some of the features. To run this demo, simply type the following
    % into the MATLAB console.
    %
    %     MarkdownPanel.demo()
    %
    % ------
    % **Attribution**
    %
    % Copyright (c) <2016> [Jonathan Suever][2].
    % All rights reserved
    %
    % This software is licensed under the [BSD license][3]
    %
    % [1]: https://github.com/showdownjs/showdown
    % [2]: https://github.com/suever
    % [3]: https://github.com/suever/MarkdownPanel/blob/master/LICENSE

    properties
        Content     = ''        % Markdown content to be displayed
        StyleSheets = {}        % List of stylesheets to link in
        Classes     = {}        % CSS classes applied to the primary div
        Options     = struct()  % Options to pass to showdown
    end

    properties (Access = 'protected')
        browser         % Handle to the javahandle_withcallbacks
        container       % Graphics handle to the HTMLBrowserPanel
        jbrowser        % Java Handle to the HTMLBrowserPanel
        listener        % Listener for when the jbrowser is deleted
        htmlComponent   % Java handle to embedded HTML component
    end

    methods
        function self = MarkdownPanel(varargin)
            % MarkdownPanel - Constructor for MarkdownPanel object
            %
            %   Can accept any property as parameter/value pair. Creates
            %   the MarkdownPanel and returns a handle to be used to
            %   manipulate the appearance/content.
            %
            % USAGE:
            %   panel = MarkdownPanel(params)
            %
            % INPUTS:
            %   params: Parameter/Value pairs, Properties to set upon
            %           creation.
            %
            % OUTPUTS:
            %   panel:  Graphics Object, The graphics handle that can be
            %           used to manipulate the appearance of the control.

            % Create an instance of the internal HTMLBrowserPanel
            import com.mathworks.mlwidgets.html.*;

            % For pre-HG2 browsers, specify the default to be HTMLPANEL
            if verLessThan('matlab', '8.4')
                HtmlComponentFactory.setDefaultType('HTMLPANEL');
            end

            % Use inputParser so we can handle structs AND params
            ip = inputParser();
            ip.KeepUnmatched = true;
            ip.addParamValue('Parent', [], @ishghandle);
            ip.parse(varargin{:});

            if ~isempty(ip.Results.Parent)
                parent = ip.Results.Parent;
            else
                parent = gcf;
            end

            % Create the Java browser component
            self.jbrowser = HTMLBrowserPanel();
            [self.browser, self.container] = javacomponent(self.jbrowser, ...
                [0 0 1 1], parent);
            self.htmlComponent = self.jbrowser.getHtmlComponent();

            % By default, make it take up the entire parent
            set(self.container, 'Units', 'norm', 'position', [0 0 1 1])

            % Now make this look like the container object by creating
            % shadow properties that interact with the underlying object
            props = fieldnames(get(self.container));

            for k = 1:numel(props)
                % Ignore if the property is already defined
                if ~isempty(self.findprop(props{k})); continue; end

                % Add a dynamic property and assign setters/getters that
                % will relay properties between the two objects
                prop = self.addprop(props{k});
                prop.SetMethod = @(s,v)setwrapper(s,prop,v);
                prop.GetMethod = @(s,e)getwrapper(s,prop);
            end

            % If the underlying graphics object is deleted, follow suit
            self.listener = addlistener(self.container, ...
                'ObjectBeingDestroyed', @(s,e)delete(self));

            % Setup the default options
            self.Options = struct(...
                'tables', true);

            % Finally consider all input arguments
            set(self, varargin{:})

            self.refresh(true);
        end

        function delete(self)
            % delete - Delete the MarkdownPanel and associated objects
            %
            % USAGE:
            %   panel.delete()

            if ishghandle(self.container)
                delete(self.container)
            end

            % Make sure that we dispose of the html component
            self.htmlComponent.dispose();
        end

        function refresh(self, force)
            % refresh - Force a refresh of the displayed markdown
            %
            % USAGE:
            %     panel.refresh(force)
            %
            % INPUTS:
            %   force:  Logical, Indicates whether to completely redraw the
            %           page (including HTML) (true) or not (false). The
            %           default is to simply execute javascript on the
            %           existing page.

            if iscell(self.Content)
                content = sprintf('%s\\n\\n', self.Content{:});
                % Remove trailing newlines
                content = regexprep(content, '\n*$', '');
            else
                content = self.Content;
            end

            % Replace "true" newlines with "\n"
            content = regexprep(content, '\n', '\\n');

            % Make sure that we properly escape double quotes so that the
            % created javascript is valid
            content = regexprep(content, '"', '\\"');

            % Javascript to run to update the HTML and make all hyperlinks
            % external
            jscript = [...
                'try {', ...
                  'var html = conv.makeHtml("', content, '");', ...
                  'display.innerHTML = html;', ...
                  'var links = document.querySelectorAll("a");', ...
                  'for (var k in links) {' ...
                    'if ( links[k].href && links[k].href.substring(0, 4) === "http" )', ...
                        ' links[k].target = "_blank";', ...
                  '}', ...
                '} catch (err) { ', ...
                  'error.innerHTML = err.message;', ...
                '}'];

            % Initial load with entire javascript
            if isempty(self.htmlComponent.getHtmlText()) || ...
                (exist('force', 'var') && force)

                % Load showdown from file that way we can catch any import
                % issues and display them in the HTML
                curdir = fileparts(mfilename('fullpath'));
                showdownjs = fullfile(curdir, 'showdown.min.js');

                % Attempt to protect the user in case they deleted showdown
                if ~exist(showdownjs, 'file')
                    % Then go download it
                    url = 'https://cdn.rawgit.com/showdownjs/showdown/1.3.0/dist/showdown.min.js';
                    urlwrite(url, showdownjs);
                end

                fid = fopen(showdownjs, 'rb');
                showdown = fread(fid, '*char');
                fclose(fid);

                % Create stylesheet entries
                if numel(self.StyleSheets)
                    format = '<link rel="stylesheet" href="%s">\n';
                    stylesheets = sprintf(format, self.StyleSheets{:});
                else
                    stylesheets = '';
                end

                % Create the options that we want to pass to the converter
                options = self.Options;
                fields = fieldnames(options);

                opts = '';

                for k = 1:numel(fields)
                    value = options.(fields{k});

                    if ischar(value)
                        value = cat(2, '"', value, '"');
                    else
                        value = num2str(value);
                    end

                    newopt = ['conv.setOption("', fields{k}, '", ', value, ');'];
                    opts = cat(2, opts, newopt);
                end

                html = {...
                    '<!DOCTYPE html>', ....
                      '<head>', ...
                        '<meta http-equiv="X-UA-Compatible" content="IE=edge">', ...
                         stylesheets, ...
                         '<script>', ...
                           'if (window.console) {', ...
                             'var console = window.console;', ...
                           '}', ...
                         '</script>', ...
                      '</head>', ...
                      '<body>', ...
                        '<div class="', sprintf('%s ', self.Classes{:}), '">', ...
                          '<div id="error" class="error" style="color:#F00"></div>', ...
                          '<div id="display">Loading...</div>', ...
                        '</div>', ...
                        '<script>', ...
                          'var display = document.getElementById("display");', ...
                          'var error = document.getElementById("error");', ...
                          'try {', ...
                            showdown(:)', ...
                            'var conv = new showdown.Converter();', ...
                            opts, ...
                          '} catch (err) {', ...
                            'error.innerHTML = err.message;', ...
                            'display.innerHTML = "";', ...
                          '}', ...
                        '</script>', ...
                      '</body>', ...
                    '</html>'};
                html = sprintf('%s\n', html{:});

                self.htmlComponent.setHtmlText(html);
            end

            % Update the HTMLBrowserPanel to use this HTML
            self.htmlComponent.executeScript(jscript);
        end
    end

    % Set/Get Methods
    methods
        function set.Content(self, val)
            % Look and see if this is a cell array
            self.Content = val;
            self.refresh();
        end

        function set.Options(self, val)
            % Check to see if they are equal to the old value
            if isequal(val, self.Options)
                return;
            end

            self.Options = val;

            % Force a complete refresh
            self.refresh(true);
        end

        function set.StyleSheets(self, val)
            if ischar(val); val = {val}; end

            self.StyleSheets = val;

            % Do a hard-refresh of the page
            self.refresh(true);
        end

        function set.Classes(self, val)
            if ischar(val); val = {val}; end

            self.Classes = val;

            % Do a hard-refresh of the page
            self.refresh(true);
        end
    end

    % These methods automatically translate the properties between the
    % underlying object and the current object
    methods (Access = 'protected')
        function setwrapper(self, prop, value)
            % Relays "set" events to the underlying container object
            set(self.container, prop.Name, value)
        end

        function res = getwrapper(self, prop)
            % Relays "get" events to the underlying container object
            res = get(self.container, prop.Name);
        end
    end

    methods (Static)
        function panel = demo()
            % demo - Demonstrate how to create/use the MarkdownPanel object
            %
            %   This demo creates a simple markdown editor/preview window
            %   that showcases how to set the content of the markdown
            %   panel. To do this, it simply shows the help text for the
            %   MarkdownPanel in both the editor and preview.
            %
            %   It also demonstrates the use of stylesheets (in this case
            %   Twitter Bootstrap)
            %
            % USAGE:
            %   panel = MarkdownPanel.demo()

            fig = figure( ...
                'Position',     [0 0 1200, 700], ...
                'Toolbar',      'none', ...
                'menubar',      'none', ...
                'NumberTitle',  'off', ...
                'Name',         'MarkdownPanel Demo');

            movegui(fig, 'center');
            drawnow;

            % Create two side-by-side panels
            flow = uiflowcontainer('v0', 'FlowDirection', 'lefttoright');

            % Grab the help section and do a little cleanup
            h = help(mfilename('fullpath'));
            h = regexprep(h, '^\s*', '');
            h = regexprep(h, '\n  ', '\n');

            % Create a java control because the builtin editbox doesn't
            % easily return the current value
            je = javax.swing.JEditorPane('text', h);
            jp = javax.swing.JScrollPane(je);
            [~, hcomp] = javacomponent(jp, [], flow);
            set(hcomp, 'Position', [0 0 0.5 1])

            % Construct the MarkdownPanel object
            twitter = 'https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap.min.css';
            panel = MarkdownPanel( ...
                'Content',      h, ...
                'Parent',       flow, ...
                'StyleSheets',  twitter, ...
                'Classes',      'container');

            % Set the option to enable smooth live previews
            panel.Options.smoothLivePreview = true;

            % Setup a timer to refresh the MarkdownPanel periodically
            timerFcn = @(s,e)set(panel, 'Content', char(je.getText()));
            htimer = timer( ...
                'Period',        1, ...
                'BusyMode',      'drop', ...
                'TimerFcn',      timerFcn, ...
                'ExecutionMode', 'fixedRate');

            % Destroy the timer when the panel is destroyed
            callback = @(s,e)delete(htimer);
            L = addlistener(panel, 'ObjectBeingDestroyed', callback);
            setappdata(fig, 'Timer', L);

            % Start the refresh timer
            start(htimer)
        end
    end
end
