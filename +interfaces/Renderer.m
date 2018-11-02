classdef Renderer <interfaces.DialogProcessor
    %RENDERER: interfaces for renderers
    % subclass has to implement renderfunction. Render is called by
    % obj.render(locs,p). locs: either localizationData or
    % localizationData.loc
    % p from Gui and inputParameters
    properties
        layer=1;
    end
    methods
        function out=run(obj,p)
             
             ax=initaxis(p.resultstabgroup,'SR image');
             image=obj.render(obj.locData,p);
             out=image;
             imagesc(image.rangex,image.rangey,image.image,'Parent',ax)
             axes(ax,'equal')

        end
        function image=render(obj,locD,p)
            if isa(locD,'interfaces.LocalizationData')
                locs=locD.getloc({'xnm','ynm', 'frame','phot','znm'},'Position','fov','layer',obj.layer);
            else
                locs=locD;
            end
            image=obj.renderfunction(locs,p);
            
        end
        function out=renderfunction(locs,p)
            out=[];
            warning('renderfunction needs to be implemented by renderer')
        end
            
    end
end
