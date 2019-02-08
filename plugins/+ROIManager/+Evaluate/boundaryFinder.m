classdef boundaryFinder<interfaces.SEEvaluationProcessor
    % This plug-in is dependent of the BALM_fibril_growth.
    % Green line is the original boundary
    % White line is the refined boundary

    properties
        boundary
    end
    methods
        function obj=boundaryFinder(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj, inp)
            out=runBoundaryFinder(obj, inp);
        end
     
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end

end



function pard=guidef(obj)
pard.t_gridRange.object=struct('Style','text','String','Bin size of the grids');
pard.t_gridRange.position=[1,1];
pard.t_gridRange.Width=2;

pard.gridRange.object=struct('Style','edit','String',5);
pard.gridRange.position=[1,3];
pard.gridRange.TooltipString = 'If you set it as 5, it means before the density comparison, every grid will be set to cover a 5-by-5 area in the original coordinates';
            
pard.t_slideStep.object=struct('Style','text','String','Slide step(s)');
pard.t_slideStep.position=[2,1];
pard.t_slideStep.Width=2;

pard.slideStep.object=struct('Style','edit','String',5);
pard.slideStep.position=[2,3];
pard.slideStep.TooltipString = 'If you set it as 5, it means during the density comparison, every grid value will be campared to its following 4 (5 minus 1, which means the reference grid itself) right neighbors';

pard.t_adjM.object=struct('Style','text','String','Adjustment of M');
pard.t_adjM.position=[3,1];
pard.t_adjM.Width=2;

pard.adjM.object=struct('Style','edit','String',1.0003);
pard.adjM.position=[3,3];
pard.adjM.TooltipString = 'If you set it as 1.0003, it means during the optimization, if the measurment of current step (Mcur) is 0.0003-time worse than the measurment of the previous step (Mpre), Mcur will still be considered as a good result. The measurment, which defines the boundary is good or not, can be definde by users.';

pard.t_mergeSteps.object = struct('Style','text','String','Order of merging steps');
pard.t_mergeSteps.position=[4,1];
pard.t_mergeSteps.Width = 2;

pard.mergeSteps.object = struct('Style','edit','String','0 3 1 2 3 1 2 0');
pard.mergeSteps.position=[4,3];
pard.mergeSteps.TooltipString = '0 = clean up; 1 = by space(ambiguous); 2 = by space(small), 3 = by time';

pard.t_minT.object = struct('Style','text','String','Minimum time');
pard.t_minT.position=[5,1];
pard.t_minT.Width = 2;

pard.minT.object = struct('Style','edit','String',5);
pard.minT.position=[5,3];
pard.minT.TooltipString = 'minimum time of a step (arbitrary unit)';

pard.t_minWidth.object = struct('Style','text','String','Minimum width');
pard.t_minWidth.position=[6,1];
pard.t_minWidth.Width = 2;

pard.minWidth.object = struct('Style','edit','String',2);
pard.minWidth.position=[6,3];
pard.minWidth.TooltipString = 'minimum width of a step (arbitrary unit)';

pard.t_method.object=struct('Style','text','String','Method (rough boundary)');
pard.t_method.position=[7,1];
pard.t_method.Width=2;

pard.method.object=struct('Style','popupmenu','String',{{'Contour line','Percentile','cumulative'}}, 'Value', 3, 'Callback', {{@mathod_callBack,obj}});
pard.method.position=[7,3];
pard.method.Width=2;

pard.t_std.object = struct('Style','text','String','Std.');
pard.t_std.position=[8,1];
pard.t_std.Width = 0.5;
pard.t_std.Visible = 'off';

pard.std.object = struct('Style','edit','String',100);
pard.std.position=[8,1.5];
pard.std.TooltipString = 'minimum width of a step (arbitrary unit)';
pard.std.Width = 0.5;
pard.std.Visible = 'off';

pard.t_contourLevel.object = struct('Style','text','String','Level');
pard.t_contourLevel.position=[8,2];
pard.t_contourLevel.Width = 0.5;
pard.t_contourLevel.Visible = 'off';

pard.contourLevel.object = struct('Style','edit','String',90);
pard.contourLevel.position=[8,2.5];
pard.contourLevel.TooltipString = 'The level (out of 100) chosen as the outline of the rough boundary';
pard.contourLevel.Width = 0.5;
pard.contourLevel.Visible = 'off';

pard.t_prctile.object = struct('Style','text','String','Prctile');
pard.t_prctile.position=[8,1];
pard.t_prctile.Width = 1;

pard.prctile.object = struct('Style','edit','String',85);
pard.prctile.position=[8,2];
pard.prctile.TooltipString = 'Percentile of dx at each time point';
pard.prctile.Width = 0.5;

% pard.dxt.Width=3;
pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi'};
pard.plugininfo.type='ROI_Evaluate';

end



function mathod_callBack(a,b,obj)
    switch obj.getSingleGuiParameter('method').Value
        case 1
            obj.guihandles.contourLevel.Visible = 'on';
            obj.guihandles.std.Visible = 'on';
            obj.guihandles.t_contourLevel.Visible = 'on';
            obj.guihandles.t_std.Visible = 'on';
            obj.guihandles.t_prctile.Visible = 'off';
            obj.guihandles.prctile.Visible = 'off';
        case 2
            obj.guihandles.contourLevel.Visible = 'off';
            obj.guihandles.std.Visible = 'off';
            obj.guihandles.t_contourLevel.Visible = 'off';
            obj.guihandles.t_std.Visible = 'off';
            obj.guihandles.t_prctile.Visible = 'on';
            obj.guihandles.prctile.Visible = 'on';
        case 3
            obj.guihandles.contourLevel.Visible = 'off';
            obj.guihandles.std.Visible = 'off';
            obj.guihandles.t_contourLevel.Visible = 'off';
            obj.guihandles.t_std.Visible = 'off';
            obj.guihandles.t_prctile.Visible = 'off';
            obj.guihandles.prctile.Visible = 'off';
    end
end