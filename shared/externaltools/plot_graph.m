function opt= plot_graph(names,from,to,varargin)
%wrapper for mwdot that plots a connection map (aka node graph)
%plot_graph(names,from,to)
%
%names : 1xn cellstring, text in each node
%from  : 1xm double,     index of start of each connection
%to    : 1xm double,     index of end   of each connection
%
%options : colour
%plot_graph(names,from,to,'-colour','rrbbkg')
%plot_graph(names,from,to,'-color' ,'rrbbkg')
%colours each box differently (use brgk)
%
%colour can also be given as an array of colours and an index:
%plot_graph(names,from,to,'-colour',[1 2 1 1],[0 0 1 ; 1 0 0])
%where the first array is the index of colours for each node
%and the second array is the colours
%
%plot_graph(names,from,to,'-comment',comments)
%adds cellstring array comments to each box
%
%plot_graph(names,from,to,'-highlight',{names})
%highlights all links from/to the given nodes with different colours
%
%plot_graph(names,from,to,'-saveas','test.pdf')
%saves to a pdf file using the most appropriate A-size paper (e.g. A4)
%saved figures are rotated and scaled to fill the page
%
%plot_graph(names,from,to,'-saveas','test.pdf','-papersize','usletter')
%saves to a pdf file using the specified paper size
%because saveas does not let you specify the fontsize that gets used, you
%may have to specify a slightly larger papersize.
%
%opionts : -shape,'rsrsrrrsr'; uses rounded/box on a per-node basis
%           default uses boxes
%
%       : -fontsize,12 %otherwise scales to fit box
%
%outputs:
%opt = plot_graph...
%returns handles to the graph created for further processing
%
%examples:
%c = {'This is A','This is B','This is C'};
%plot_graph({'A','B','C'},[1 1 3],[1 2 1],'-colour','rgb','-comment',c);
%
%acknowledge : mGraphViz (file exchange ID 24652)
%acknowledge : mkdotfile (file exchange ID 27608)
%
%uses the local file plot_graph.mwdot to temporarily store figure data
%this is overwritten with every use.

opt = sub_opt(names,from,to,varargin{:});
opt = sub_mdot_write(opt);
opt = sub_mwdot_call(opt);
opt = sub_figure_create(opt);
opt = sub_node_create(opt);
opt = sub_comment_add(opt);
opt = sub_fontsize(opt);
opt = sub_edge_create(opt);
opt = sub_arrow_create(opt);
opt = sub_saveas(opt);
end

function opt = sub_opt(names,from,to,varargin)
opt.names = names;
for i=1:numel(names)
    str = ['A' num2str(i)];
    while numel(str)<numel(names{i}); str = ['A' str]; end
    opt.namemw{i} = str;
end
opt.deps.from = from;
opt.deps.to   = to;
opt.node_col = repmat('b',size(names));
opt.edge_col = repmat('k',size(from));
opt.node_round = repmat('r',size(names));
opt.comments = {};
opt.fig.axes = [];
opt.fig.main = [];
opt.saveas = '';
opt.papersize = '';
opt.highlight = {};
opt.fontsize = 'auto';
opt.direction= 'UD';
i = 0;
while i<numel(varargin)
    i = i+1;
    val = varargin{i};
    switch val
        case {'-colours','-colors','-colour','-color'};
            opt.node_col = varargin{i+1}; i = i+1;
            if isnumeric(opt.node_col)
                opt.node_col_array = varargin{i+1}; i = i+1;
            end
        case {'-shape'}; opt.node_round = varargin{i+1}; i=i+1;
        case {'-comments','-comment'}; opt.comments = varargin{i+1}; i=i+1;
        case {'-print','-saveas'}; opt.saveas = varargin{i+1}; i=i+1;
        case {'-papersize'}; opt.papersize = varargin{i+1}; i = i+1;
        case {'-highlight'}; opt.highlight = varargin{i+1}; i = i+1;
        case {'-fontsize'}; opt.fontsize = varargin{i+1}; i = i+1;
        case {'-edgeColor'};
            opt.edge_col=varargin{i+1}; i = i+1;
        case {'-direction'}
            opt.direction=varargin{i+1}; i = i+1;
        otherwise
            if ischar(val); opt.col_node = val; continue; end
            if iscellstr(val); opt.comments = val; continue; end
            if ishghandle(val); opt.fig.axes = val; continue; end
            disp(val); error('unknown option');
    end
end
%check sizes
if ~isequal(length(opt.names),length(opt.node_col));
    error('colours are wrong size');
end
if ~isequal(size(opt.names),size(opt.comments)) && numel(opt.comments)
    error('comments are wrong size');
end

if numel(opt.node_round)==1; opt.node_round = repmat(opt.node_round,size(names)); end

%process colours to array, if they aren't
if ischar(opt.node_col)
    opt.node_col_array = unique(opt.node_col)';
    [~,opt.node_col] = ismember(opt.node_col,opt.node_col_array);
end

%name of file to save the list of node connections to
%this will be the input to mwdot
opt.mwdot.infile = [which(mfilename) 'dot'];
opt.mwdot.infile = strrep(opt.mwdot.infile,'\','/');
if exist(opt.mwdot.infile,'file'); delete(opt.mwdot.infile); end
end

function opt = sub_mdot_write(opt)
%create commands for mdot file
str{1} = sprintf('digraph %s_calls{',opt.namemw{1});
direc=upper(opt.direction);
str{2} = ['rankdir="' direc '"'];
% str{2} = 'rankdir="LR"';
str{3} = 'node [shape=box, color=blue, facecolor=white]';
for i=1:numel(opt.deps.from)
    str{end+1} = sprintf('%s -> %s;',...
        opt.namemw{opt.deps.from(i)},...
        opt.namemw{opt.deps.to(i)});
end
str{end+1} = '}';

%write mdot file
fid = fopen(opt.mwdot.infile,'wt');
for i=1:numel(str)
    fprintf(fid,'%s\n',str{i});
end
fclose(fid);
end

function opt = sub_mwdot_call(opt)
%calls the mwdot function that ships with matlab
%feeds it the mdot file we prepared earlier


%added 20150309 need som additional logic here to catch all variants of os
%if you do get a problem here, please get in touch so I can add your configuration

fol = lower(computer('arch'));
switch fol;
    case 'win32';   foo = 'mwdot.exe';  %not tested
    case 'win64';   foo = 'mwdot.exe';  %tested
    case 'maci64';  foo = 'mwdot';      %not tested
    case 'glnxa64'; foo = 'mwdot';      %not tested
    case 'glnx86';  foo = 'mwdot';      %not tested
    otherwise;
        %if you get this error, please get in touch so I can add your OS !
        error('do not recognise this OS');
end

opt.mwdot.file = fullfile(matlabroot,'bin',fol,foo);

%if you get this error, please get in touch
if ~exist(opt.mwdot.file,'file')
    error('did not find file "%s"',strrep(opt.mwdot.file,'\','/'));
end

%call executable
%20150309 added " quotes around filenames, to handle space in foldernames
opt.mwdot.call = ['"' opt.mwdot.file '" -Tplain "' opt.mwdot.infile '"'];
[status,result]=system(opt.mwdot.call);

%throw an error if no engine was found
if status~=0
    disp(opt.mwdot.call);
    error([foo ' failed with error "' result '"']);
end
opt.mwdot.result = result;
end

function opt = sub_figure_create(opt)
%opt = sub_figure_create(opt)
%creates the figure and axis

opt.fig.main = gcf;
clf;
opt.fig.axes = gca;

%parse the "graph" line
pat = ['^graph\s+' ,...
    '(?<size>[^\s]+)\s+' ,...
    '(?<wide>[^\s]+)\s+' ,...
    '(?<high>[^\s]+)\s+'];
val = regexp(opt.mwdot.result,pat,'tokens','once');
num = str2num(char(val)); %#ok<ST2NM>
opt.fig.scale = num(1);
opt.fig.shape = 1;   %x/y shape ratio
opt.fig.wide  = num(2);
opt.fig.high  = num(3);

%adjust the scale so it will fit in the available screen size
inchSize = get(0,'ScreenSize')/get(0,'ScreenPixelsPerInch');
inchSize = inchSize([3 4]);
screen.wide = inchSize(1)-1; %leave a 1/2" margin left/right
screen.high = inchSize(2)-1.5; %3/4" margin top/bottom
opt.fig.scale = min([screen.wide/opt.fig.wide screen.high/opt.fig.high 1]);

% opt.fig.scale=1;
%axes fill figure
set(opt.fig.axes,...
    'units','normalized',...
    'position',[0 0 1 1],...
    'XTick',[],...
    'YTick',[],...
    'XLim',opt.fig.wide*opt.fig.scale*opt.fig.shape*[-0.01 1.01],...
    'YLim',opt.fig.high*opt.fig.scale*[-0.01 1.01])
hold on;

%scale figure so axes are right shape
set(opt.fig.main,...
    'units','inches',...
    'position',[0 0 opt.fig.wide*opt.fig.scale*opt.fig.shape opt.fig.high*opt.fig.scale],...
    'units','pixels');
centerfig(opt.fig.main,0);
axis image

%set the figure bacground and axes colors to axes background
set(gcf,'Color' ,get(gca,'Color'));
set(gca,'YColor',get(gca,'Color'))
set(gca,'XColor',get(gca,'Color'))
end

function opt = sub_node_create(opt)

%this pattern straight from mGraphViz
pat = ['node\s+' ...
    '(?<nodeID>"[^"]+"|[\w\d]+)\s+' ...
    '(?<x>[\d\.\+e]+)\s+' ...
    '(?<y>[\d\.\+e]+)\s+' ...
    '(?<w>[\d\.\+e]+)\s+' ...
    '(?<h>[\d\.\+e]+)\s+' ...
    '(?<label>".*"|[\w\d]+)'];
nodes = regexp(opt.mwdot.result,pat,'names');
%nodes.x = str2num(char({nodes.x}));
%TODO use rectangle command instead of patch
for i=1:numel(nodes)
    opt.node(i).text =            nodes(i).label;
    opt.node(i).x    = str2double(nodes(i).x)*opt.fig.scale*opt.fig.shape;
    opt.node(i).y    = str2double(nodes(i).y)*opt.fig.scale;
    opt.node(i).w    = str2double(nodes(i).w)*opt.fig.scale*opt.fig.shape;
    opt.node(i).h    = str2double(nodes(i).h)*opt.fig.scale;
    [tf,loc] = ismember(opt.node(i).text,opt.namemw);
    if ~tf; keyboard; end
    opt.node(i).text  = opt.names{loc};
    opt.node(i).index = loc;
    opt.node(i).color = opt.node_col_array(opt.node_col(loc),:);
    opt.node(i).round = opt.node_round(loc);
    
    %create with 10 point fontsize as standard
    %this may get adjusted in sub_fontsize
    opt.node(i).handle_text = text(...
        opt.node(i).x,...
        opt.node(i).y,...
        opt.node(i).text,...
        'HorizontalAlignment','center',...
        'Interpreter','none','verticalAlignment','middle',...
        'margin',20,...
        'edgeColor','none',...
        'FontSize',10);
    
        opt.node(i).handle_box = rectangle(...
            'Position',...
            [opt.node(i).x-opt.node(i).w/2 opt.node(i).y-opt.node(i).h/2 opt.node(i).w opt.node(i).h],...
            'EdgeColor',opt.node(i).color);
        switch opt.node(i).round
            case 's'; c = [0 0];
            case 'r'; c = [opt.node(i).h/opt.node(i).w 1];
            case 'e'; c = [1 1];
        end
        set(opt.node(i).handle_box,'Curvature',c);
        
end
opt = rmfield(opt,'node_col');
opt = rmfield(opt,'node_col_array');
opt = rmfield(opt,'node_round');

end

function opt = sub_fontsize(opt)
%fit text to boxes with 10% margin

fscale = inf;
if ~isnumeric(opt.fontsize)
    for i=1:numel(opt.node)
        pos_txt = get(opt.node(i).handle_text,'Extent');
        pos_box = get(opt.node(i).handle_box,'Position');
        fscale = min([fscale 0.9*pos_box(3)/pos_txt(3) 0.9*pos_box(4)/pos_txt(4)]);
        
    end
    %round down to nearest .5, do not go below 6 point
    fs = fscale*get(opt.node(1).handle_text,'FontSize');
    fs = max([6 fs]);
    opt.fontsize = floor(fs*2)/2;
end

for i=1:numel(opt.node)
    set(opt.node(i).handle_text,'FontSize',opt.fontsize);
end
end

function opt = sub_edge_create(opt)
%opt = sub_edge_create(opt)
%creates lines connecting nodes

pat = ['edge\s+' ...
    '(?<tail>(?:"[^"]+")|(?:[\w\d]+))\s+' ...
    '(?<tip>(?:"[^"]+")|([\w\d]+))\s+' ...
    '(?<numSpline>\d+)\s+'...
    '(?<coords>[\d\.\+e\s]+)\s+'];
edges = regexp(opt.mwdot.result,pat,'names');
for i=1:numel(edges)
    xy = str2num(edges(i).coords);
    x = xy(1:2:end)*opt.fig.shape*opt.fig.scale;
    y = xy(2:2:end)*opt.fig.scale;
    opt.edge(i).endpoints = {edges(i).tail edges(i).tip};
    [tf,loc] = ismember(opt.edge(i).endpoints,opt.namemw);
    if ~all(tf); keyboard; end
    opt.edge(i).index = loc;
    opt.edge(i).endpoints = opt.names(opt.edge(i).index);
    
    %correct for lines that go right to left
    try
        %which nodes does this point to ?
        [tf,loc] = ismember(opt.edge(i).index,[opt.node.index]);
        if ~all(tf); keyboard; end
        if diff([opt.node(loc).x])<0;
            x = x(end:-1:1);
            y = y(end:-1:1);
        end
    catch
        keyboard
    end
    %these points are bezier guide points, not points on the actual line
    %make a bezier curve. Note this curve does not go through the points
    [~,~,val] = sub_bezier(0,1,[x ; y]',0:.01:1);
    opt.edge(i).handle = plot(val(:,1),val(:,2),opt.edge_col(i));
end
opt = rmfield(opt,'edge_col');

%if highlight option, find all links to/from highlights nodes and colour them
val = [opt.edge.index];
[tf,loc] = ismember(opt.highlight,opt.names);
loc = loc(tf);
if numel(loc)
    for i = 1:numel(opt.edge)
        set(opt.edge(i).handle,'Color',[.85 .85 .85]);
    end
    for i = find(ismember(val(2:2:end),loc))
        set(opt.edge(i).handle,'Color',[0 .5 1]);
    end
    for i = find(ismember(val(1:2:end),loc))
        set(opt.edge(i).handle,'Color',[1 0 0]);
    end
end
end

function opt = sub_arrow_create(opt)
%opt = sub_arrow_create(opt)
%creates arrows at the end of edges

%scale arrows by height of text
scale = get(opt.node(1).handle_text,'extent');
%if there are more than one lines of text, divide as appropriate
str = get(opt.node(1).handle_text,'String');
if iscell(str); scale = scale/numel(str); end

for i=1:numel(opt.edge)
    %angle of arrow from last pair of points on line
    x = get(opt.edge(i).handle,'XData'); %x = x((end-1):end);
    y = get(opt.edge(i).handle,'YData');% y = y((end-1):end);
    [~,ind]=min(y);
    if strcmpi(opt.direction,'ud')
        if ind==1
    x=x([2 1]);y=y([2 1]);
        else
            x=x(end-1:end);y=y(end-1:end);
        end
    else
        x=x(end-1:end);y=y(end-1:end);
    end
    ang = atan2(diff(y),diff(x))+pi/2;
    
    xo = .4*[0  1  0 -1]*scale(4);
    yo = .4*[0 .6 -1 .6]*scale(4);
%    plot(xo,yo);
    
    xs = xo*cos(ang) - yo*sin(ang);
    ys = xo*sin(ang) + yo*cos(ang);
    
    col = get(opt.edge(i).handle,'Color');
    opt.edge(i).handle_arrow = patch(...
        x(end)+xs,...
        y(end)+ys,...
        zeros(size(xs)),...
        'edgecolor',col,...
        'facecolor',col);
end

end %function sub_arrow_create

function opt = sub_comment_add(opt)
if isempty(opt.comments); return; end
for i = 1:numel(opt.node)
    h = opt.node(i).handle_text;
    set(h,'String',{get(h,'String'),opt.comments{i}})
end
end

function opt = sub_saveas(opt)
if isempty(opt.saveas); opt = rmfield(opt,{'saveas','papersize'}); return; end



orient(opt.fig.main,'portrait');
%20140908 scale figure to fit paper, not screen
%szf = opt.fig.scale*[opt.fig.wide*opt.fig.shape opt.fig.high]; %size of figure
szf = [opt.fig.wide opt.fig.high];

if isempty(opt.papersize)
    szp = [0 0]; A = 6;
    while szp(1)<min(szf) || szp(2)<max(szf)
        A = A-1;
        set(opt.fig.main,'PaperType',['A' num2str(A)]);
        szp = get(opt.fig.main,'PaperSize');
        if A==0; break; end
    end
    opt.papersize = get(opt.fig.main,'PaperType');
else
    set(opt.fig.main,'PaperType',opt.papersize);
end

if szf(1)>szf(2); orient(gcf,'landscape');
else              orient(gcf,'portrait');
end

%scale x and y to fit page
szp = get(opt.fig.main,'PaperSize');
scale = szf./szp;
scale = max(scale);
%place figure in center of paper
offx = szp(1)/2-szf(1)/2/scale;
offy = szp(2)/2-szf(2)/2/scale;
set(opt.fig.main,'PaperPosition',[offx offy szf(1)/scale szf(2)/scale]);

%turn the filename into a proper filename with full path
[fol,foo,ext] = fileparts(opt.saveas);
if isempty(fol); fol = pwd; end
if isempty(ext); ext = '.pdf'; end
if ~isequal(lower(ext),'.pdf'); error('can only save as pdf'); end
opt.saveas = fullfile(fol,[foo ext]);
saveas(opt.fig.main,opt.saveas);
end

function [X,Y,val] = sub_bezier(a,b,p,y)
% function val = bezier(a,b,p,y)
% bezier(0,1,,linspace(0,1,nplot))
% INPUT: a = left edge
%        b = right edge,
%        p = control points
%        y = series of linear points to use for bezier.
%
% OUTPUT: Bezier polynomial
% Courtesy of Jonas Ballani

n = size(p,1);
m = length(y);
val = zeros(m,2);
X(:,1) = p(:,1);
Y(:,1) = p(:,2);

for j = 1:m
    for i = 2:n
        X(i:n,i) = (b-y(j))/(b-a)*X(i-1:n-1,i-1) + (y(j)-a)/(b-a)*X(i:n,i-1);
        Y(i:n,i) = (b-y(j))/(b-a)*Y(i-1:n-1,i-1) + (y(j)-a)/(b-a)*Y(i:n,i-1);
    end
    val(j,1) = X(n,n);
    val(j,2) = Y(n,n);
end
end %function bezier

function sub_not_used()
end

%% DEVNOTES
%branched from mGraphViz because:
%(1) it does not work under 2014B or mac
%(2) I can get rid of a lot of mwdot options that I don't use
%(3) can control the size of the containing figure
%(4) allows be to keep it self-contained with plot_depfun

%done 20140904 scale figure to fit screen and use equal size axis
%done 20140904 accept input to automatically save to pdf
%done 20140904 remove "axes" option, force use of entire current figure
%done 20140905 highlight all edges to/from a given node
%done 20140904 saving to figure sort of works, but can be unreliable.
%              limited to pdf for now, as this seems most robust
%done 20140905 bezier function as subroutine (cleared with original author)
%done 20150121 margin to prevent top of figure falling off top of screen
%current version 20150121
%todo additional colour "h" for "hidden" is gray.
%todo option to make boxes bigger with less space between them
%todo force given fontsize
