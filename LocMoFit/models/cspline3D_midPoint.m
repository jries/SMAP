classdef cspline3D_midPoint<geometricModel
    % A c-spline for describing a linear structure traversing in 3D space.
    %
    % Geometric parameters:
    %   * `xMid, yMid, zMid`: (nm) the xyz coordinates of the mid point.
    %   * `dist`: (nm): the distance between neighbouring control points
    %   * `rotAzi[L/R]_n, rotEle[L/R]_n`: (°) the azimuthal and elevation
    % angles or of the vector pointing to the [L/R]_n control
    % point. `[L/R]` is either L (left) or R (right) with respect to the
    % mid point. `n` indicates the order. For exmple, rotAziL1 means the
    % azimuthal angle defining the 1st point on the left of the mid
    % point.
    %
    % Relavent biological structure:
    %   * Actin filaments
    %   * The central axes of microtubules
    %
    % Preview:
    %   .. image:: ./images/models/cspline3D_midPoint.PNG
    %       :width: 400
    %   Scale bar: 50 nm.
    methods
        function obj = cspline3D_midPoint(varargin)
            obj@geometricModel(varargin{:});
            % Define the default argument values here in the constructor.
            obj.name = {'xMid', 'yMid', 'zMid', 'dist','rotAziL1','rotEleL1', 'rotAziR1','rotEleR1','rotAziL2','rotEleL2', 'rotAziR2','rotEleR2'};
            obj.fix = [0 0 0 0 0 0 0 0 0 0 0 0] ;
            obj.value = [0 0 0 50 0 0 0 0 0 0 0 0];
            obj.lb = [-30 -30 -30 0 -45 -45 -45 -45 -45 -45 -45 -45];
            obj.ub = [30 30 30 200 45 45 45 45 45 45 45 45];
            obj.min = [-500 -500 -500 -inf -inf -inf -inf -inf -inf -inf -inf -inf];
            obj.max = [500 500 500 inf inf inf inf inf inf inf inf inf];


            obj.modelType = 'discretized';
            obj.modelTypeOption = {'discretized', 'continuous'};
            obj.dimension = 3;
            obj.internalSettings.numOfCtrlPointSet = 2;
            obj.listed = true;
        end
        function [model, p]= reference(obj, par, dx)
        % For details, see :meth:`reference`.

            % set additional parameters of the model
            ip = inputParser;
            fn = fieldnames(obj.internalSettings);
            for k = 1:length(fn)
                ip.addParameter(fn{k}, obj.internalSettings.(fn{k}));
            end
            a = {};
            parse(ip,a{:})
            pResults = ip.Results;
            
            numOfCtrlPointSet = pResults.numOfCtrlPointSet;
            
            % control points:
            [ctrlX,ctrlY,ctrlZ] = getCtrlPoints(par, numOfCtrlPointSet);
            
            if isempty(obj.ParentObject.locsPrecFactor)
                locsPrecFactor = 1;
            else
                locsPrecFactor = obj.ParentObject.locsPrecFactor;
            end
            
            minD = locsPrecFactor*dx;
            
            % vecterization
            for l = 1:size(ctrlZ,2)
                arcLen = arclength(ctrlX(:,l), ctrlY(:,l), ctrlZ(:,l),'linear');
                samplingFactor = round(arcLen/minD);
                [pt,diffVector] = interparc(double(samplingFactor), double(ctrlX(:,l)), double(ctrlY(:,l)), double(ctrlZ(:,l))) ;
            end
            
            model.x = pt(:,1);
            model.y = pt(:,2);
            model.z = pt(:,3);
            model.channel = ones(size(model.x));
            model.n = ones(size(model.x));
            
            %     likelihood = gaussDist(model,double(locs.ynmrot), double(locs.xnmrot), double(locs.znm), 8,8,10);
            %     likelihoodBoun = gaussDist(model,0,0,0, 16,16,20);
            
            p.numOfCtrlPointSet=numOfCtrlPointSet;
        end
        
        function items = getThings2Plot(obj,par)
        % The user can define what should be also displayed in the plots.
        % Args:
        % 	obj: a :class:`geometricModel` object.
        % Returns:
        %   items: things to be plotted.
        %
            numOfCtrlPointSet = obj.getInternalSettings('numOfCtrlPointSet');
            [ctrlX,ctrlY,ctrlZ] = getCtrlPoints(par, numOfCtrlPointSet);
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
                lRotAziL = startsWith(obj.name,'rotAziL');
                if any(lRotAziL)
                    ctrlSetID = str2double(regexprep(obj.name(lRotAziL),'\D',''));
                    lastCtrlPoint_ori = max(ctrlSetID);
                else
                    lastCtrlPoint_ori = 0;
                end
                obj.ParentObject.ParentObject.resetInit;
                % detect the addition or removement of control points
                if lastCtrlPoint_ori<obj.getInternalSettings('numOfCtrlPointSet')
                    % when the last control point number is smaller than the
                    % numOfCtrlPoint
                    for k = lastCtrlPoint_ori+1:obj.getInternalSettings('numOfCtrlPointSet')
                        if k ~= 0||lastCtrlPoint_ori~=obj.getInternalSettings('numOfCtrlPointSet')
                            obj.addMPar(defaultParsArgNewCtrlPoint(k));
                        end
                    end
                elseif obj.getInternalSettings('numOfCtrlPointSet')<lastCtrlPoint_ori
                    % when the last control point number is larger than the
                    % numOfCtrlPoint
                    obj.ParentObject.ParentObject.resetInit;
                    for k = obj.getInternalSettings('numOfCtrlPointSet')+1:lastCtrlPoint_ori
                        if lastCtrlPoint_ori~=obj.getInternalSettings('numOfCtrlPointSet')
                            obj.rmMPar(defaultParsArgNewCtrlPoint(k).name);
                        end
                    end
                end
                % notify the SMLMModel obj to update the internal parameters
                notify(obj.ParentObject,'mParsArgModified');
                evtdata = mParsArgModifiedData(obj.ParentObject.ID);
                notify(obj.ParentObject.ParentObject,'mParsArgModified',evtdata);
                if ~isempty(obj.ParentObject.ParentObject.linkedGUI)
                    notify(obj.ParentObject.ParentObject.linkedGUI,'mParsArgModified',evtdata);
                end
                obj.ParentObject.ParentObject.saveInit;
            end
        end
    end
end

function struct = defaultParsArgNewCtrlPoint(number)
struct.name = {['rotAziL' num2str(number)],['rotEleL' num2str(number)],['rotAziR' num2str(number)],['rotEleR' num2str(number)]};
struct.fix = [0 0 0 0] ;
struct.value = [0 0 0 0];
struct.lb = [-12.5 -12.5 -12.5 -12.5];
struct.ub = [12.5 12.5 12.5 12.5];
struct.min = [-inf -inf -inf -inf];
struct.max = [inf inf inf inf];
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

function [ctrlX,ctrlY,ctrlZ] = getCtrlPoints(par, numOfCtrlPointSet)
% control points:
    xMid = par.xMid;
    yMid = par.yMid;
    zMid = par.zMid;
    dist = par.dist;
    
    % left side
    for k = 1:numOfCtrlPointSet
        rotAziL(k,:) = par.(['rotAziL' num2str(k)]); % ->
        rotEleL(k,:) = par.(['rotEleL' num2str(k)]); % ->
        rotAziR(k,:) = par.(['rotAziR' num2str(k)]); % ->
        rotEleR(k,:) = par.(['rotEleR' num2str(k)]); % ->
    end
    rotAziL = cumsum(rotAziL,1);
    rotEleL = cumsum(rotEleL,1);
    rotAziL = deg2rad(rotAziL);
    rotEleL = deg2rad(rotEleL);
    [diffXL,diffYL,diffZL] = sph2cart(rotAziL,rotEleL,-dist);
    
    rotAziR = cumsum(rotAziR,1);
    rotEleR = cumsum(rotEleR,1);
    rotAziR = deg2rad(rotAziR);
    rotEleR = deg2rad(rotEleR);
    [diffXR,diffYR,diffZR] = sph2cart(rotAziR,rotEleR,dist);
    
    %% 191128
    diffXL = cumsum([xMid; diffXL],1);
    diffYL = cumsum([yMid; diffYL],1);
    diffZL = cumsum([zMid; diffZL],1);
    
    diffXR = cumsum([xMid; diffXR],1);
    diffYR = cumsum([yMid; diffYR],1);
    diffZR = cumsum([zMid; diffZR],1);
    
    ctrlX = [diffXL(end:-1:2);diffXR];
    ctrlY = [diffYL(end:-1:2);diffYR];
    ctrlZ = [diffZL(end:-1:2);diffZR];
end