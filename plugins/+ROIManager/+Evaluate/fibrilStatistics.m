classdef fibrilStatistics<interfaces.SEEvaluationProcessor
    % This is a plugin in development. Public has no access to the 
    % run_ functions called in this plugin. For internal users, 
    % "fibrilKymograph" is required.
    methods
        function obj=fibrilStatistics(varargin)        
            obj@interfaces.SEEvaluationProcessor(varargin{:});
            flagDirExist = exist('../fibrilkymograph','dir');
            if flagDirExist==0
                warning('Please install fibrilkymograph.')
            else
                addpath(genpath('../fibrilkymograph'))
            end
        end
        
        function out=run(obj,p)
            if p.filtering
                if isfield(obj.site.evaluation,'fibrilStatistics')&&isfield(obj.site.evaluation.fibrilStatistics,'setting')
                else
                    p.axisLb = 0;
                    p.axisUb = 0;
                end
            end
            out = runFibrilStatistics(obj,p);
            
            if isfield(obj.site.evaluation,'fibrilStatistics')&&isfield(obj.site.evaluation.fibrilStatistics,'setting')
                obj.guihandles.axisLb.String = out.setting.axisLb;
                obj.guihandles.axisUb.String = out.setting.axisUb;
            else
                obj.guihandles.axisLb.String = 0;
                obj.guihandles.axisUb.String = obj.site.evaluation.fibrilStraightener.straightenedSA.arclen;
            end
            if ~p.lockFilter
                obj.guihandles.filtering.Value = 0;
                obj.guihandles.msg.Visible = 0;
            end
            
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function pard = guidef(obj)



    pard.viewOnly.object = struct('Style','checkbox','String','View only', 'Value', 0, 'Callback', {{@viewOnlyCallback,obj}});
    pard.viewOnly.position = [1 1];
    pard.viewOnly.Width = 1;
    
    pard.t_edgeFactor.object = struct('Style','text','String','Edge: +-');
    pard.t_edgeFactor.position = [1 3];
    pard.t_edgeFactor.Width = 1;
    
    pard.edgeFactor.object = struct('Style','edit','String','100');
    pard.edgeFactor.position = [1 4];
    pard.edgeFactor.Width = 1;
    
    pard.t1.object = struct('Style','text','String','Filtering:');
    pard.t1.position = [2 1];
    pard.t1.Width = 1;
    
    pard.axisLb.object = struct('Style','edit','String','0');
    pard.axisLb.position = [2 2];
    pard.axisLb.Width = 1;
    
    pard.axisUb.object = struct('Style','edit','String','0');
    pard.axisUb.position = [2 3];
    pard.axisUb.Width = 1;
    
    pard.filtering.object = struct('Style','checkbox','String','','Value',0, 'Callback', {{@filterCallback,obj}});
    pard.filtering.position = [2 4];
    pard.filtering.Width = 0.3;
    
    pard.lockFilter.object = struct('Style','checkbox','String','Lock','Value',0, 'Callback', {{@lockCallback,obj}});
    pard.lockFilter.position = [2 4.3];
    pard.lockFilter.Width = 1;
    pard.lockFilter.Enable = 'off';

    pard.resetLength.object = struct('Style','pushbutton','String','Reset length','Value', 0, 'Callback', {{@resetLenCallback,obj}});
    pard.resetLength.position = [3 1];
    pard.resetLength.Width = 1;
    
    pard.msg.object = struct('Style','text','String','Click on the site to save change!');
    pard.msg.position = [3 2];
    pard.msg.Width = 3;
    pard.msg.Visible = 'off';
    
    pard.tFrame.object = struct('Style','text','String','Frame:');
    pard.tFrame.position = [4 1];
    pard.tFrame.Width = 1;
    
    pard.frameLb.object = struct('Style','edit','String','0');
    pard.frameLb.position = [4 2];
    pard.frameLb.Width = 1;
    
    pard.frameUb.object = struct('Style','edit','String',max(obj.locData.loc.frame));
    pard.frameUb.position = [4 3];
    pard.frameUb.Width = 1;
    
    pard.old.object = struct('Style','checkbox','String','Old version','Value', 0);
    pard.old.position = [5 1];
    pard.old.Width = 1;
    
    pard.locsSource.object = struct('Style','popupmenu','String',{{'Smooth_axis','Rough_axis'}},'Value', 1);
    pard.locsSource.position = [5 2];
    pard.locsSource.Width = 1;

    % check whether these is a fibrilDynamics before this module
    options = obj.locData.SE.processors.eval.guihandles.modules.Data(:,2);      % names of all mounted modules
    lFibrilDynamics = contains(options,'fibrilDynamics');                   % all fibrilDynamics modules
    if sum(lFibrilDynamics)==0
        pard.frameLb.Enable = 'off';
        pard.frameUb.Enable = 'off';
    end    
    
    pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi'};
    pard.plugininfo.type='ROI_Evaluate';
end

function filterCallback(a,b,obj)
    if obj.guihandles.filtering.Value
        obj.guihandles.lockFilter.Enable = 'on';
        obj.guihandles.msg.Visible = 'on';
    else
        obj.guihandles.lockFilter.Enable = 'off';
        obj.guihandles.lockFilter.Value = 0;
        obj.guihandles.msg.Visible = 'off';
    end
end

function resetLenCallback(a,b,obj)
    if isfield(obj.locData.SE.currentsite.evaluation, 'fibrilStraightener')
        obj.guihandles.axisUb.String = obj.locData.SE.currentsite.evaluation.fibrilStraightener.straightenedSA.arclen;
        obj.guihandles.axisLb.String = 0;
        obj.guihandles.filtering.Value = 1;
        obj.guihandles.msg.Visible = 'on';
        obj.guihandles.lockFilter.Enable = 'off';
        obj.guihandles.lockFilter.Value = 0;
    else
        disp('Please run "fibrilStraightener" first')
    end
end

function lockCallback(a,b,obj)
    if obj.guihandles.lockFilter.Value
        obj.guihandles.msg.Visible = 'off';
    elseif obj.guihandles.filtering.Value
        obj.guihandles.msg.Visible = 'on';
    else
        obj.guihandles.msg.Visible = 'off';
    end
end

function viewOnlyCallback(a,b,obj)
    fn = fieldnames(obj.guihandles);
    if a.Value == 1
        for k = 1:length(fn)
            obj.guihandles.(fn{k}).Visible = 'off';
        end
    else
        for k = 1:length(fn)
            obj.guihandles.(fn{k}).Visible = 'on';
        end
    end
    obj.guihandles.viewOnly.Visible = 'on';
    lockCallback([],[],obj)
end
