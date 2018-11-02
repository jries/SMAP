function mymsgbox(varargin)
h=msgbox(varargin{:});
h.Visible='off';
txth= findobj(h,'Type','text','-depth',3);
xo=txth.Extent;
txth.FontSize=16;
xn=txth.Extent;
h.Position=h.Position+xn-xo;
h.Visible='on';