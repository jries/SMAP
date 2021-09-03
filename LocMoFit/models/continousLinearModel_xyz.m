classdef continousLinearModel_xyz<geometricModel
    % A general model describing a linear structure traverse in 3D space.
    methods
        function obj = continousLinearModel_xyz(varargin)
            obj@geometricModel(varargin{:});
            % Define the default argument values here in the constructor.
            obj.name = {'x0', 'y0', 'z0', 'x1','y1','z1', 'x2', 'y2','z2', 'x3', 'y3','z3'};
            obj.fix = [0 0 0 0 0 0 0 0 0 0 0 0] ;
            obj.value = [0 0 0 100 0 0 200 0 0 300 0 0];
            obj.lb = [-30 -30 -30 -30 -150 -150 -50 -150 -150 -50 -150 -150];
            obj.ub = [30 30 30 50 150 150 50 150 150 30 150 150];
            obj.min = [-500 -500 -500 -500 -500 -500 -500 -500 -500 -500 -500 -500];
            obj.max = [500 500 500 500 500 500 500 500 500 500 500 500];
            obj.modelType = 'discretized';
            obj.modelTypeOption = {'discretized', 'continuous'};
            obj.dimension = 3;
            obj.internalSettings.numOfCtrlPoint = 4;
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
struct.name = {['x' num2str(number)], ['y' num2str(number)],['z' num2str(number)]};
struct.fix = [0 0 0] ;
struct.value = [0 0 0];
struct.lb = [-100 -100 -100];
struct.ub = [100 100 100];
struct.min = [-500 -500 -500];
struct.max = [500 500 500];
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

ctrlX = zeros(numOfCtrlPoint,1);
ctrlY = zeros(numOfCtrlPoint,1);
ctrlZ = zeros(numOfCtrlPoint,1);

ctrlX(1,1) = x0;
ctrlY(1,1) = y0;
ctrlZ(1,1) = z0;

for k = 1:numOfCtrlPoint-1
    ctrlX(k+1,:) = par.(['x' num2str(k)]); % ->
    ctrlY(k+1,:) = par.(['y' num2str(k)]); % ->
    ctrlZ(k+1,:) = par.(['z' num2str(k)]); % ->
end
end