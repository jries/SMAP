classdef continousLinearModel<geometricModel
    % A general model describing a linear structure traverse in 3D space.
    methods
        function obj = continousLinearModel(varargin)
            obj@geometricModel(varargin{:});
            % Define the default argument values here in the constructor.
            obj.name = {'x0', 'y0', 'z0', 'dist1','rotAzi1','rotEle1', 'dist2', 'rotAzi2','rotEle2', 'dist3', 'rotAzi3','rotEle3'};
            obj.fix = [0 0 0 0 0 0 0 0 0 0 0 0] ;
            obj.value = [0 0 0 50 0 0 50 0 0 50 0 0];
            obj.lb = [-30 -30 -30 0 -45 -45 0 -45 -45 0 -45 -45];
            obj.ub = [30 30 30 200 45 45 200 45 45 200 45 45];
            obj.min = [-500 -500 -500 1 -inf -inf 1 -inf -inf 1 -inf -inf];
            obj.max = [500 500 500 400 inf inf 400 inf inf 400 inf inf];
            obj.modelType = 'discretized';
            obj.modelTypeOption = {'discretized', 'continuous'};
            obj.dimension = 3;
            obj.internalSettings.numOfCtrlPoint = 3;
        end
        function [model, p]= reference(obj, par, dx)
            % set additional parameters of the model
            ip = inputParser;
            fn = fieldnames(obj.internalSettings);
            for k = 1:length(fn)
                ip.addParameter(fn{k}, obj.internalSettings.(fn{k}));
            end
            a = {};
            parse(ip,a{:})
            pResults = ip.Results;
            
            numOfCtrlPoint = pResults.numOfCtrlPoint;
            
            % control points:
            [ctrlX,ctrlY,ctrlZ] = getCtrlPoints(par, numOfCtrlPoint);
            samplingFactor = round(dx);
            
            % vecterization
            for l = 1:size(ctrlZ,2)
                [pt,diffVector] = interparc(20, ctrlX(:,l), ctrlY(:,l), ctrlZ(:,l)) ;
            end
            
            model.x = pt(:,1);
            model.y = pt(:,2);
            model.z = pt(:,3);
            model.channel = ones(size(model.x));
            model.n = ones(size(model.x));
            
            %     likelihood = gaussDist(model,double(locs.ynmrot), double(locs.xnmrot), double(locs.znm), 8,8,10);
            %     likelihoodBoun = gaussDist(model,0,0,0, 16,16,20);
            
            p.numOfCtrlPoint=numOfCtrlPoint;
        end
        
        function items = getThings2Plot(obj,par)
        % The user can define what should be also displayed in the plots.
        % Args:
        % 	obj: a :class:`geometricModel` object.
        % Returns:
        %   items: things to be plotted.
        %
            numOfCtrlPoint = obj.getInternalSettings('numOfCtrlPoint');
            [ctrlX,ctrlY,ctrlZ] = getCtrlPoints(par, numOfCtrlPoint);
            items(1).XData = ctrlX;
            items(1).YData = ctrlY;
            items(1).ZData = ctrlZ;
            items(1).Marker = 'x';
            items(1).MarkerSize = 12;
            items(1).MarkerEdgeColor = 'b';
            items(1).LineStyle = 'none';
            items(1).LineWidth = 3;
        end
    end
    methods(Access = protected)
        function respond2InternalSettingsChange(obj,event)
            if isempty(obj.ParentObject)||isempty(obj.ParentObject.ParentObject)
            else
                lastCtrlPoint_ori = str2double(regexprep(obj.name{end},'\D',''));
                % detect the addition or removement of control points
                if lastCtrlPoint_ori<obj.getInternalSettings('numOfCtrlPoint')-1
                    % when the last control point number is smaller than the
                    % numOfCtrlPoint-1
                    obj.ParentObject.ParentObject.resetInit;
                    for k = lastCtrlPoint_ori+1:obj.getInternalSettings('numOfCtrlPoint')-1
                        if k ~= 0||lastCtrlPoint_ori~=obj.getInternalSettings('numOfCtrlPoint')-1
                            obj.addMPar(defaultParsArgNewCtrlPoint(k));
                        end
                    end
                elseif obj.getInternalSettings('numOfCtrlPoint')-1<lastCtrlPoint_ori
                    % when the last control point number is larger than the
                    % numOfCtrlPoint-1
                    obj.ParentObject.ParentObject.resetInit;
                    for k = obj.getInternalSettings('numOfCtrlPoint'):lastCtrlPoint_ori+1
                        if lastCtrlPoint_ori~=obj.getInternalSettings('numOfCtrlPoint')-1
                            obj.rmMPar(defaultParsArgNewCtrlPoint(k).name);
                        end
                    end
                end
                % notify the SMLMModel obj to update the internal parameters
                notify(obj.ParentObject,'mParsArgModified');
                evtdata = mParsArgModifiedData(obj.ParentObject.ID);
                notify(obj.ParentObject.ParentObject,'mParsArgModified',evtdata);
                notify(obj.ParentObject.ParentObject.linkedGUI,'mParsArgModified',evtdata);
                obj.ParentObject.ParentObject.saveInit;
            end
        end
    end
end

function struct = defaultParsArgNewCtrlPoint(number)
struct.name = {['dist' num2str(number)], ['rotAzi' num2str(number)],['rotEle' num2str(number)]};
struct.fix = [0 0 0] ;
struct.value = [50 0 0];
struct.lb = [0 -45 -45];
struct.ub = [200 45 45];
struct.min = [1 -inf -inf];
struct.max = [400 inf inf];
end

function T = getRotMat(p1,t)
p0 = [1 0 0];
p1 = p1./sqrt(p1*p1');
C = cross(p0', p1') ;
D = dot(p0', p1') ;
NP0 = norm(p0) ; % used for scaling
if ~all(C==0) % check for colinearity
    Z = [0 -C(3) C(2); C(3) 0 -C(1); -C(2) C(1) 0] ;
    R = (eye(3) + Z + Z^2 * (1-D)/(norm(C)^2)) / NP0^2 ; % rotation matrix
else
    R = sign(D) * (norm(p1) / NP0) ; % orientation and scaling
end
T = [[R;[0 0 0]] [t'; 1]];
end

function [ctrlX,ctrlY,ctrlZ] = getCtrlPoints(par, numOfCtrlPoint)
% control points:
    x0 = par.x0;
    y0 = par.y0;
    z0 = par.z0;

    for k = 1:numOfCtrlPoint-1
        rotAzi(k,:) = par.(['rotAzi' num2str(k)]); % ->
        rotEle(k,:) = par.(['rotEle' num2str(k)]); % ->
        dist(k,:) = par.(['dist' num2str(k)]); % ->
    end
    rotAzi = cumsum(rotAzi,1);
    rotEle = cumsum(rotEle,1);
    rotAzi = rotAzi.*pi./180;
    rotEle = rotEle.*pi./180;

    [diffX,diffY,diffZ] = sph2cart(rotAzi,rotEle,dist);
    %% 191128
    diffX = [x0;diffX];
    diffY = [y0;diffY];
    diffZ = [z0;diffZ];
    ctrlX = cumsum(diffX,1);
    ctrlY = cumsum(diffY,1);
    ctrlZ = cumsum(diffZ,1);
end