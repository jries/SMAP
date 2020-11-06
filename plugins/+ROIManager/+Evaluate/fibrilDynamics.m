classdef fibrilDynamics<interfaces.SEEvaluationProcessor
    % This plug-in is dependent of the BALM_fibril_growth.
    % Green line is the original boundary
    % White line is the refined boundary
    
    % log
    % 
    
    properties
        poly
        hpoly
        boundary
    end
    methods
        function obj=fibrilDynamics(varargin)
            obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj, inp)
            out=runFibrilDynamics(obj, inp);
            valh=obj.site.evaluation.(obj.modulename);
            if length(obj.poly)>=obj.site.ID &&~isempty(obj.poly{obj.site.ID})
                out.manualBound = obj.dynamicsManualBound;
            elseif isfield(valh,'poly')
                obj.poly{obj.site.ID}=valh.poly;
                out.manualBound = obj.dynamicsManualBound;
            end
        end
        
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        
        function addline(obj,a,b,h)
            % if length(obj.poly)<obj.site.ID || isempty(obj.poly{obj.site.ID})
            obj.hpoly=impoly(h,'Closed',false,'PositionConstraintFcn',@obj.polyconstrain);
            % else %add vortex to end
            %     pos=obj.hpoly.getPosition;
            %     pos(end+1,:)=pos(end,:);
            %     obj.hpoly.setPosition(vertcat(pos(1,:),pos));
            
            %  end
            setConstrainedPosition(obj.hpoly,obj.hpoly.getPosition);
            hi=addNewPositionCallback(obj.hpoly,@obj.polycallback);
            obj.polycallback(obj.hpoly.getPosition);

        end
        
        function rmline(obj,a,b,h)
            % 200929: this funciton is created as an option for removing
            % the drawn line
            answer = questdlg('Remove the line?');
            switch answer
                case 'Yes'
                case 'No'
                    return
                case 'Cancel'
                    return
            end
            child1 = h.Children(1);
            if isa(child1,'matlab.graphics.primitive.Group')
                delete(child1)
                obj.poly{obj.site.ID}=[];
                obj.site.evaluation.(obj.name) = rmfield(obj.site.evaluation.(obj.name), 'poly');
            else
                disp('Noting to remove.')
            end
        end
        
        function outp=polyconstrain(obj,inp)
            
            outp=vertcat(inp(1,:),max(inp(2:end,:),inp(1:end-1,:)));
            
        end
        
        function polycallback(obj,inp)
%             setConstrainedPosition(obj.hpoly,inp);
            obj.poly{obj.site.ID}=inp;
            obj.site.evaluation.(obj.name).poly=inp;
            obj.site.evaluation.fibrilDynamics.manualBound = obj.dynamicsManualBound;
        end
        
        function output = dynamicsManualBound(obj)
            % This function gets the steps from segments
            
            %% Basic info.
            % these should match the same parameters defined in the histogram.
            binFactorPos = 10;
            binFactorTime = 2.5;
            
            % this value is newly defined here:
            stallThreshold = 0.1;
            
            %% Get segements
            % manualBound.data: anchor points of the segmentation
            % manualBound.segment: segments between anchor points
            % manualBound.ratePerSeg: diffPos/diffTime per segment
            
            manualBound.data = obj.poly{obj.site.ID}; % column: pos time
            manualBound.data = [manualBound.data(:,1)*binFactorPos manualBound.data(:,2)*binFactorTime];
            manualBound.segment = diff(manualBound.data); % segments that form the boundary
            manualBound.ratePerSeg = manualBound.segment(:,1)./manualBound.segment(:,2);
            
            %% Robust measures
            lStall = manualBound.ratePerSeg<stallThreshold;
            fullPosNTime = manualBound.data(end,:)-manualBound.data(1,:);
            % total growth time
            growthTime = sum(manualBound.segment(~lStall,2));
            % total growth length
            growthLength = sum(manualBound.segment(~lStall,1));
            % total stall time/total time
            stallFraction = 1-(growthTime/fullPosNTime(2));
            % total growth length/total growth time
            growthOnlyRate = growthLength/growthTime;
            
            % Get steps. Multiple segments in a row will be merged into one
            % step.
            % 'mark' denotes growth/stall phases:
            %   *positive: growth
            %   *negative: stall
            % Individual phases have their own IDs.
            %
            mark = cumsum(lStall)+1;
            mark(lStall)=-mark(lStall);
            
            manualBound.ratePerSeg = [manualBound.ratePerSeg mark];
            
            if mark(end) < 0
                % remove the last step if it is a stall
                manualBound.segment(end,:) = [];
                manualBound.ratePerSeg(end,:) = [];
                mark(end,:) = [];
            end
            
            orderSteps = unique(mark,'stable');       % This is the order of the steps

%             avgRate = sum(manualBound.data,1);
            
            avgRate = fullPosNTime(1)/fullPosNTime(2);
            
            [sumPos,names]=grpstats(manualBound.segment(:,1),mark,{'sum','gname'});
            [sumTime,~]=grpstats(manualBound.segment(:,2),mark,{'sum','gname'});
            stepMark = str2num(char(names));
            timePerStall = sumTime(lStall);
            
            [~,ind] = ismember(orderSteps, stepMark);
            manualBound.steps = [sumPos(ind) sumTime(ind)];
            firstTimePoint = manualBound.data(1,2);
            
            endtimeStep = firstTimePoint + cumsum(manualBound.steps(:,2));
            midTimeStep = endtimeStep - manualBound.steps(:,2)/2;
            manualBound.steps = [manualBound.steps midTimeStep orderSteps];
            
%             validation = firstTimePoint+sum(manualBound.steps(:,2)); %    this value should be the same end time of the manual boundary.
%             validation
            
            if stepMark(1) < 0
                % if the first grp is the stall, then rm it from the
                % sumTime and sumPos.
                lStall = stepMark<0;
                stallTime = sum(sumTime(lStall));
                sumTime(lStall) = [];
                sumPos(lStall) = [];
            else
                stallTime = 0;
            end
            stepRate = sumPos./sumTime;
            
            %% Export
            % values per fibril
            output.steps = manualBound.steps;
            output.stallTime = stallTime;
            output.avgRate = avgRate;
            output.growthTime = growthTime;
            output.stallFraction = stallFraction;
            output.growthOnlyRate = growthOnlyRate;
            
            % values per growth phase
            output.stepRate = stepRate;
            output.stepWidth = sumPos;
            output.stepSpan = sumTime;

            % values per stall phase
            output.timePerStall = timePerStall;
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

pard.t_DistFScale.object=struct('Style','text','String','DistF scale');
pard.t_DistFScale.position=[3,1];
pard.t_DistFScale.Width=1;

pard.DistFScale.object=struct('Style','edit','String',7e-2);
pard.DistFScale.position=[3,2];
pard.DistFScale.TooltipString = 'Lower this value to tolerate slower increase of foreground density (Df)';
pard.DistFScale.Width=1;

pard.t_offSetScale.object=struct('Style','text','String','Offset scale');
pard.t_offSetScale.position=[3,3];
pard.t_offSetScale.Width=1;

pard.offSetScale.object=struct('Style','edit','String',1e-1);
pard.offSetScale.position=[3,4];
pard.offSetScale.TooltipString = 'Higher this value to allow higher background density';

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

pard.displayOpt.object = struct('Style','checkbox','String','show optimisation');
pard.displayOpt.position=[7,1];
pard.displayOpt.Width = 2;

pard.isFlip.object = struct('Style','checkbox','String','Flip(position)');
pard.isFlip.position=[7,3];
pard.isFlip.Width = 1.5;

pard.isFlipSave.object = struct('Style','checkbox','Callback', {{@isFlipSave_callback,obj}});
pard.isFlipSave.position=[7,4.5];
pard.isFlipSave.Width = 0.2;

pard.isFlipLock.object = struct('Style','checkbox');
pard.isFlipLock.position=[7,4.7];
pard.isFlipLock.Width = 0.2;
pard.isFlipLock.Enable = 'off';

pard.t_frameCut.object = struct('Style','text', 'string', 'Frame cut:');
pard.t_frameCut.position=[8,1];
pard.t_frameCut.Width = 1;

pard.frameCut.object = struct('Style','edit', 'string', 0);
pard.frameCut.position=[8,2];
pard.frameCut.Width = 1;

pard.t_frameCut.object = struct('Style','text', 'string', 'Frame cut:');
pard.t_frameCut.position=[8,1];
pard.t_frameCut.Width = 1;

pard.isManual.object = struct('Style','checkbox','string', 'Manual', 'Value',1);
pard.isManual.position=[9,1];
pard.isManual.Width = 1;

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

function isFlipSave_callback(a,b,obj)
if ~obj.guihandles.isFlipSave.Value
    obj.guihandles.isFlipLock.Value = 0;
    obj.guihandles.isFlipLock.Enable = 'off';
else
    obj.guihandles.isFlipLock.Enable = 'on';
end
end