classdef csplineClosedTube3D_midPoint<geometricModel
    % A c-spline for describing a flexible tube that has its ends closed, traversing in 3D space.
    %
    % Geometric parameters:
    %   * `xMid, yMid, zMid`: (nm) the xyz coordinates of the mid point.
    %   * `r`: (nm) the radius of the tube.
    %   * `dist`: (nm): the distance between neighbouring control points
    %   * `rotAzi[L/R]_n, rotEle[L/R]_n`: (°) the azimuthal and elevation angles or of the vector pointing to the [L/R]_n control point. `[L/R]` is either L (left) or R (right) with respect to the mid point. `n` indicates the order. For exmple, rotAziL1 means the azimuthal angle defining the 1\ :sup:`st` point on the left of the mid point.
    %
    % Relavent biological structure:
    %   * outer membrane of a mitochondria
    %
    % Preview:
	% 	.. note::
	% 		It will be available soon.
	%
	% ..
    %   .. image:: ./images/models/csplineClosedTube3D_midPoint.PNG
    %       :width: 400
    %   Scale bar: 50 nm.
    methods
        function obj = csplineClosedTube3D_midPoint(varargin)
            obj@geometricModel(varargin{:});
            % Define the default argument values here in the constructor.
            obj.name = {'xMid', 'yMid', 'zMid','r', 'dist','rotAziL1','rotEleL1', 'rotAziR1','rotEleR1','rotAziL2','rotEleL2', 'rotAziR2','rotEleR2'};
            obj.fix = [0 0 0 0 0 0 0 0 0 0 0 0 0] ;
            obj.value = [0 0 0 30 50 0 0 0 0 0 0 0 0];
            obj.lb = [-30 -30 -30 -inf -inf -15 -15 -15 -15 -15 -15 -15 -15];
            obj.ub = [30 30 30 inf inf 15 15 15 15 15 15 15 15];
            obj.min = [-500 -500 -500 10 1 -inf -inf -inf -inf -inf -inf -inf -inf];
            obj.max = [500 500 500 60 400 inf inf inf inf inf inf inf inf];

            obj.modelType = 'discretized';
            obj.modelTypeOption = {'discretized', 'continuous'};
            obj.dimension = 3;
            obj.internalSettings.numOfCtrlPointSet = 2;
            obj.listed = true;
        end
        function [model, p]= reference(obj, par, dx)
        % For details, see :meth:`reference`.

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
            
            for l = 1:size(ctrlZ,2)% vecterization
                arcLen = arclength(ctrlX(:,l), ctrlY(:,l), ctrlZ(:,l),'linear');
                samplingFactor = round(arcLen/minD);
                [pt,diffVector] = interparc(samplingFactor, ctrlX(:,l), ctrlY(:,l), ctrlZ(:,l)) ;
            end
            
            if par.r > 0
                % put rings perpenticular to the axis of the curve
                nOnRing = round(2*pi*par.r/minD);
                theta = 0:2*pi/nOnRing:(2*pi-2*pi/nOnRing);
                % build a ring based on the radius par.r
                [x,y] = pol2cart(theta,par.r);
                xyRing = [zeros(size(x')) x' y'];
                tube = zeros(nOnRing, 3, [size(diffVector,1)]-2);
                point_hs = [];
                for k = 1:size(diffVector,1)
                    if k == size(diffVector,1)
                        T = getRotMat(diffVector(k,:),pt(k,:));
                    else
                        T = getRotMat(-diffVector(k,:),pt(k,:));
                    end
                    if ismember(k, [1 size(diffVector,1)])
                        % add hemisphere to the ends of the tube
                        [x,y,z] = hemisphere(par.r, locsPrecFactor,dx);
                        temp = T*[z; y; x; ones(size(z))];
                        point_hs= [point_hs; temp(1:end-1,:)'];
                    else
                        temp = T*[xyRing'; ones(1,length(xyRing))];
                        tube(:,:,k-1)= temp(1:end-1,:)';
                    end
                end

                tube = permute(tube, [1 3 2]);
                tube = reshape(tube, [],3);
                tube = [tube; point_hs];
                
                model.x = tube(:,1);
                model.y = tube(:,2);
                model.z = tube(:,3);
            else
                model.x = pt(:,1);
                model.y = pt(:,2);
                model.z = pt(:,3);
            end
            model.channel = ones(size(model.x));
            model.n = ones(size(model.x));
            
            %     likelihood = gaussDist(model,double(locs.ynmrot), double(locs.xnmrot), double(locs.znm), 8,8,10);
            %     likelihoodBoun = gaussDist(model,0,0,0, 16,16,20);
            
            p.numOfCtrlPointSet=numOfCtrlPointSet;
        end
        
        function items = getThings2Plot(obj,par, varargin)
        % The user can define what should be also displayed in the plots.
		%
        % Args:
        % 	obj: a :class:`@geometricModel.geometricModel` object.
        %
		% Returns:
        %   items: things to be plotted.
            numOfCtrlPointSet = obj.getInternalSettings('numOfCtrlPointSet');
            [ctrlX,ctrlY,ctrlZ] = getCtrlPoints(par, numOfCtrlPointSet);
            items(1).XData = ctrlX;
            items(1).YData = ctrlY;
            items(1).ZData = ctrlZ;
            items(1).Marker = 'x';
            items(1).MarkerSize = 12;
            items(1).MarkerEdgeColor = 'b';
            items(1).Color = 'none';
            items(1).LineStyle = 'none';
            items(1).LineWidth = 3;
            
            if isempty(obj.ParentObject.locsPrecFactor)
                locsPrecFactor = 1;
            else
                locsPrecFactor = obj.ParentObject.locsPrecFactor;
            end
            minD = locsPrecFactor*0.75;
            arcLen = arclength(ctrlX, ctrlY, ctrlZ,'linear');
            samplingFactor = round(arcLen/minD);
            pt = interparc(round(samplingFactor/2), ctrlX, ctrlY, ctrlZ);
            
            items(2).XData = pt(:,1);
            items(2).YData = pt(:,2);
            items(2).ZData = pt(:,3);
            items(2).Marker = 'none';
            items(2).MarkerSize = 1;
            items(2).MarkerEdgeColor = 'w';
            items(2).Color = 'w';
            items(2).LineStyle = '-';
            items(2).LineWidth = 2;
        end
        
        function derivedPars = getDerivedPars(obj, pars)
            % Exports a empty variable when no derived parameters.
            par = obj.ParentObject.ParentObject.exportPars(1,'mPar');
            % control points:
            [ctrlX,ctrlY,ctrlZ] = getCtrlPoints(par, obj.internalSettings.numOfCtrlPointSet);
            
            if isempty(obj.ParentObject.locsPrecFactor)
                locsPrecFactor = 1;
            else
                locsPrecFactor = obj.ParentObject.locsPrecFactor;
            end
            
            minD = locsPrecFactor*0.75;
            
            arcLen = arclength(ctrlX, ctrlY, ctrlZ,'spline');
            samplingFactor = round(arcLen/minD);
            pt = interparc(round(samplingFactor/2), ctrlX, ctrlY, ctrlZ);
            
            [L,R,K] = curvature(pt);
            derivedPars.curvature = 1./R;
            derivedPars.arcLength = arcLen;
            derivedPars.curvatureVector = K;
            
            derivedPars.avgCurvature = mean(derivedPars.curvature,'omitnan');
            derivedPars.p75Curvature = prctile(derivedPars.curvature, 75);
%             figure; plot3(pt(:,1),pt(:,2),pt(:,3))
%             % plot3([pt(:,1) pt(:,1)+50000*K(:,1)]',[pt(:,2) pt(:,2)+50000*K(:,2)]',[pt(:,3) pt(:,3)+50000*K(:,3)]')
%             axis equal
%             hold on; arrow3(pt,pt+K.*50000,'-k',0.2,0.4)


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
                notify(obj.ParentObject.ParentObject.linkedGUI,'mParsArgModified',evtdata);
                obj.ParentObject.ParentObject.saveInit;
            end
        end
    end
end

function struct = defaultParsArgNewCtrlPoint(number)
struct.name = {['rotAziL' num2str(number)],['rotEleL' num2str(number)],['rotAziR' num2str(number)],['rotEleR' num2str(number)]};
struct.fix = [0 0 0 0] ;
struct.value = [0 0 0 0];
struct.lb = [-inf -inf -inf -inf];
struct.ub = [inf inf inf inf];
struct.min = [-12.5 -12.5 -12.5 -12.5];
struct.max = [12.5 12.5 12.5 12.5];
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
    %     R = sign(D) * (norm(p1) / NP0) ; % orientation and scaling
    R = sign(D).*eye(3);
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

function [x,y,z] = hemisphere(r, locsPrecFactor, dx)
closingAngle = 0;
samplingFactor = (33.9/(locsPrecFactor*dx))^2;      % change this if you need more points
samplingScale = round((4*pi*r^2)/(4*pi)*0.01);
if samplingScale == 0
    samplingScale = 1;
end
samples = 1 * samplingFactor * samplingScale;       % number of sampled points.

% Unit Fibonacci sphere
rnd = 1;
offset = 2./samples;
gr = (sqrt(5.0) + 1.0) / 2.0;                       % golden ratio = 1.6180339887498948482
ga = (2.0 - gr) * (2.0*pi);                         % golden angle = 2.39996322972865332
increment = ga;
samplesUb = round((sin(closingAngle)+1)./offset-0.5);
k = 1:samplesUb;

if abs(r) < inf
    z = (k+0.5).* offset - 1;
    rs = zeros(size(k));
    if samplesUb>samples/2
        rs(1:(round(samples/2))) = (1 - z(1:(round(samples/2))).^2).^0.5;
        rs(round(samples/2)+1:samplesUb) = rs(round(samples/2):-1:2*round(samples/2)-samplesUb+1);
    else
        rs(1:samplesUb) = (1 - z(1:samplesUb).^2).^0.5;          % reduced the computational cost by
    end
    phi = rem((k + rnd), samples) .* increment; %
    
    x = cos(phi) .* rs(k);
    y = sin(phi) .* rs(k);
    
    y = real(y);
end

%% Get the final model
% Scale the unit sphere to radius
y = real(r*y);
x = real(r*x);
z = -r*z;      % upside-down the coordinates
end