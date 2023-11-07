classdef fibrilStraightener<interfaces.SEEvaluationProcessor
    % This plug-in is dependent of the BALM_fibril_growth.
    % Green line is the original boundary
    % White line is the refined boundary

    properties
        boundary
    end
    methods
        function obj=fibrilStraightener(varargin)        
            obj@interfaces.SEEvaluationProcessor(varargin{:});
            flagDirExist = exist('../fibrilkymograph','dir');
            if flagDirExist==0
                warning('Please install fibrilkymograph.')
            else
                addpath(genpath('../fibrilkymograph'))
            end
        end
        function out=run(obj, inp)
            out=runFibrilStraightener(obj, inp);
        end
     
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end

end



function pard=guidef(obj)
    % pard.dxt.Width=3;
    % pard.isFlip.object = struct('Style','checkbox','String','Flip(position)');

    pard.flip.object = struct('Style','pushbutton','String', 'Flip','Callback', {{@flip_callback,obj}});
    pard.flip.position=[1,1];
    
    pard.t_smoothFactor.object = struct('Style','text','String', 'Smooth factor');
    pard.t_smoothFactor.position=[2,1];

    pard.smoothFactor.object = struct('Style','edit','String', 1e-5);
    pard.smoothFactor.position=[2,2];

    pard.t_window.object = struct('Style','text','String', 'Window');
    pard.t_window.position=[3,1];

    pard.window.object = struct('Style','edit','String', 500);
    pard.window.position=[3,2];

    pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi'};
    pard.plugininfo.type='ROI_Evaluate';

end

function flip_callback(a,b,obj)
    if isempty(obj.site)
        warning('Not evaluated. You have to click on one site first.')
        return
    end
    if ~isfield(obj.site.evaluation.(obj.modulename), 'flip')
        disp('Flipping applied.')
        obj.site.evaluation.(obj.modulename).flip = true;
    else
        obj.site.evaluation.(obj.modulename).flip = ~obj.site.evaluation.(obj.modulename).flip;
        if obj.site.evaluation.(obj.modulename).flip
            disp('Flipping applied.')
        else
            disp('Flipping canceled.')
        end
    end
end