function guiStyle(guihandles, fn, varargin)
    % This function converts relative positions in the unit of lines
    % into absolute positions.
    p = inputParser;
    p.addParameter('exclude',[])
    p.addParameter('mode','default')
    p.addParameter('fontsize',10)
    p.addParameter('FieldHeight',26)
    p.addParameter('tabsize1',[0    -1  546 342])
    p.addParameter('tabsize2',[0    -1  540 311])
    p.addParameter('Vsep',3)
    p.addParameter('Xrim',3)
    p.addParameter('Vrim',2)
    p.parse(varargin{:});
    p = p.Results;
    switch p.mode
        case 'default'
            if ispc
                fontsize=9;
                FieldHeight=26;
                Vsep=3;
                Xrim=7;
                Vrim=15;
            elseif ismac
                fontsize=15;
                FieldHeight=25;
                Vsep=1;
                Xrim=2;
                Vrim=0;
            elseif isunix
                fontsize=9;
                FieldHeight=26;
                Vsep=3;
                Xrim=7;
                Vrim=2;
            end
        case 'userDefined'
            fontsize=9;
            FieldHeight=26;
            Vsep=3;
            Xrim=7;
            Vrim=2;
    end
	FieldWidth = FieldHeight*4;
    Xsep = 1;
    
    lKept = ~strcmp(p.exclude, fn);
    fn = fn(lKept);
    for k = 1:length(fn)
        h = guihandles.(fn{k});
        originalUnits_hP = h.Parent.Units;
        h.Parent.Units = 'pixels';
        
        edgeLen_hP = h.Parent.Position(3:4);
        pos = h.Position;
        
%         pos(1) = Xrim + (pos(1)-1)*(Xsep + FieldWidth);
        pos(1) = Xrim + (pos(1)-1)*(Xsep + FieldWidth);
%         pos(2) = edgeLen_hP(2) - (Vrim + (pos(2)-pos(4)+1)*(Vsep+FieldHeight));
        pos(2) = edgeLen_hP(2) - (Vrim + (pos(2)+pos(4))*(Vsep+FieldHeight));
        pos(3:4) = [pos(3)*FieldWidth + floor(pos(3)-1)*Xsep pos(4)*FieldHeight + floor(pos(4)-1)*Vsep];
        
        h.Units = 'pixels';
        h.FontSize = fontsize;
        h.Position = pos;
        h.Parent.Units = originalUnits_hP;
        if isa(h,'matlab.ui.control.UIControl')&&strcmp(h.Style,'text')
            h.HorizontalAlignment = 'left';
        end
    end
end