classdef NPC3Dfitting_forVilma<interfaces.SEEvaluationProcessor
    % This is a plugin in development. Public has no access to the 
    % run_ functions called in this plugin. For internal user, you need
    % NPC3D to run this plugin.
    methods
        function obj=NPC3Dfitting_forVilma(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        
        function out=run(obj,p)
            out = runNPC3DfittingJ_forVilma(obj,p);
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function pard = guidef(obj)
    pard.showFig.object = struct('Style','checkbox','String','Show fig', 'Value', 0);
    pard.showFig.position = [1 1];
    pard.showFig.Width = 1;
    
    pard.t1.object = struct('Style','text','String', 'Bounds (xpos ypos zpos distTwoRings rotXY rotXZ):');
    pard.t1.position = [2 1];
    pard.t1.Width = 4;
    
    pard.t2.object = struct('Style','text','String', 'Lower:');
    pard.t2.position = [3 1];
    pard.t2.Width = 0.7;
    
    pard.lb.object = struct('Style','edit','String', '-110 -110 -250 -100 0 -30');
    pard.lb.position = [3 1.7];
    pard.lb.Width = 3;
    
    pard.t3.object = struct('Style','text','String', 'Upper:');
    pard.t3.position = [4 1];
    pard.t3.Width = 0.7;
    
    pard.ub.object = struct('Style','edit','String', '110 110 300 0 360 30');
    pard.ub.position = [4 1.7];
    pard.ub.Width = 3;
    
    pard.t4.object = struct('Style','text','String', 'Offset:');
    pard.t4.position = [5 1];
    pard.t4.Width = 0.7;
    
    pard.sumOffset.object = struct('Style','edit','String', '1e-90');
    pard.sumOffset.position = [5 1.7];
    pard.sumOffset.Width = 1;
    
    pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi'};
    pard.plugininfo.type='ROI_Evaluate';
end