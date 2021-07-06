function h = errorshade(varargin)
    nargs = length(varargin);
    if nargs < 3
        error(message('Not enough inputs'))
    end
    if isa(varargin{1},'matlab.graphics.axis.Axes')
    	ax = varargin{1};
        varargin(1) = [];
    end
    x = varargin{1};
    y = varargin{2};
    error = varargin{3};
    if length(varargin) > 3
        p = inputParser;
        p.addParameter('Color',[]);
        p.addParameter('LineWidth',[]);
        p.addParameter('Marker',[]);
        p.addParameter('MarkerFaceColor',[]);
        p.addParameter('MarkerSize',3);
%         p.addParameter('Shade_FaceColor',[]);
        p.addParameter('Shade_EdgeColor','none');
        p.addParameter('Shade_FaceAlpha',0.5);
        p.addParameter('Shade_LineWidth',[]);
        p.addParameter('Shade_FaceColor',[]);
        p.addParameter('Shade_LineStyle',[]);
        p.parse(varargin{4:end});
        p = p.Results;
    end
    h = plot(ax,x,y);
    
    if length(varargin) > 3
        fn = fieldnames(p);
        if length(varargin) > 3
            for k = 1:length(fn)
                oneF = p.(fn{k});
                if ~isempty(oneF)&&~startsWith(fn{k}, 'Shade_')
                    h.(fn{k}) = p.(fn{k});
                end
            end
        end
    end
    hold(ax,'on')
    ub = y+error;
    lb = y-error;
    errorY = [ub;lb(end:-1:1)];
    errorX = [x;x(end:-1:1)];
    if strcmp(p.Shade_FaceColor,'none')
        pa_1 = plot(ax,x,ub);
        pa_2 = plot(ax,x,lb);
        if length(varargin) > 3
            for k = 1:length(fn)
                oneF = p.(fn{k});
                oneFnOld = fn{k};
                if ~isempty(oneF)&&startsWith(oneFnOld, 'Shade_')&&~any(strcmp(oneFnOld,{'Shade_FaceAlpha','Shade_FaceColor'}))
                    oneFnNew = replace(oneFnOld,'Shade_','');
                    if strcmp(oneFnNew, 'EdgeColor')
                        oneFnNew = 'Color';
                    end
                    pa_1.(oneFnNew) = p.(oneFnOld);
                    pa_2.(oneFnNew) = p.(oneFnOld);
                end
            end
        end
    else
        pa = fill(ax,errorX,errorY,'r');
        if length(varargin) > 3
            for k = 1:length(fn)
                oneF = p.(fn{k});
                oneFnOld = fn{k};
                if ~isempty(oneF)&&startsWith(oneFnOld, 'Shade_')
                    oneFnNew = replace(oneFnOld,'Shade_','');
                    pa.(oneFnNew) = p.(oneFnOld);
                end
            end
        end
        %     hold(ax,'off')
        if isempty(p.Shade_FaceColor)
            pa.FaceColor = p.Color;
        end
    end
    
end