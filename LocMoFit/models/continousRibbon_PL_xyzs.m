classdef continousRibbon_PL_xyzs<geometricModel
    % A general model describing a linear structure traverse in 3D space.
    methods
        function obj = continousRibbon_PL_xyzs(varargin)
            obj@geometricModel(varargin{:});
            % Define the default argument values here in the constructor.
            obj.name = {'r','cx1','cy1','cz1','s1', 'cx2', 'cy2','cz2','s2', 'cx3', 'cy3','cz3','s3'};
            obj.fix = [1 0 0 0 0 0 0 0 0 0 0 0 0] ;
            obj.value = [65 0 0 0 0 100 0 0 0 200 0 0 0];
            obj.lb = [-10 -50 -50 -50 -pi -50 -50 -50 -pi -50 -50 -50 -pi];
            obj.ub = [20 50 50 50 pi 50 50 50 50 pi 50 50 50 pi];
            obj.min = [30 -5000 -5000 -5000 -pi -5000 -5000 -5000 -pi -5000 -5000 -5000 -pi];
            obj.max = [150 5000 5000 5000 inf 5000 5000 5000 inf 5000 5000 5000 inf];
            obj.modelType = 'discretized';
            obj.modelTypeOption = {'discretized', 'continuous'};
            obj.dimension = 3;
            obj.internalSettings.numOfCtrlPointSet = 3;
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
            
            numOfCtrlPointSet = pResults.numOfCtrlPointSet;
            
            % control points:
            [ctrlX,ctrlY,ctrlZ,ctrlS] = getCtrlPointsXYZS(par, numOfCtrlPointSet);

            if isempty(obj.ParentObject.locsPrecFactor)
                locsPrecFactor = 1;
            else
                locsPrecFactor = obj.ParentObject.locsPrecFactor;
            end
            
            minD = locsPrecFactor*1; %%0.75 dx
            
            for l = 1:size(ctrlZ,2)% vecterization
                arcLen = arclength(ctrlX(:,l), ctrlY(:,l), ctrlZ(:,l),'linear');
                samplingFactor = round(arcLen/minD);
                [pt,diffVector] = interparc(samplingFactor, ctrlX(:,l), ctrlY(:,l), ctrlZ(:,l),'pchip') ; %samplingFactor
            end
            
            if par.r > 0
                % put axes perpenticular to the central line of the
                % SC/ribbon
                nOnRing = 2;
                theta = [0 pi];
                axes = zeros(nOnRing, 3, [size(diffVector,1)]);
                %%FIX:Do notconnect the angle between segment - keep  them
                %%indep
                distances=pdist2([ctrlX ctrlY ctrlZ],pt,'euclidean');
                [local_minD,local_idx] = min(distances,[],2);
                ptS=zeros(1,size(diffVector,1));
                for i_seg=1:(size(local_idx,1)-1)
                    ptS(1,local_idx(i_seg):local_idx(i_seg+1))=ctrlS(i_seg);
                end
                %ptS=spline(1:size(ctrlS,1),ctrlS,linspace(1,size(ctrlS,1),size(diffVector,1)));
                for k = 1:size(diffVector,1)
                    [x,y] = pol2cart(theta+(ptS(k)),par.r);%pi*/180
                    xyRing = [zeros(size(x')) x' y'];
                    T = getTransmat(diffVector(k,:),pt(k,:));
                    temp = T*[xyRing.'; ones(1,size(xyRing,1))];
                    axes(:,:,k)= temp(1:end-1,:)';
                end
                axes = permute(axes, [1 3 2]);
                axes = reshape(axes, [],3);
                model.x = axes(:,1);
                model.y = axes(:,2);
                model.z = axes(:,3);
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
        % Args:
        % 	obj: a :class:`geometricModel` object.
        % Returns:
        %   items: things to be plotted.
        %
            numOfCtrlPointSet = obj.getInternalSettings('numOfCtrlPointSet');
            [ctrlX,ctrlY,ctrlZ,ctrlS] = getCtrlPointsXYZS(par, numOfCtrlPointSet);
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
            minD = locsPrecFactor*1;%%0.75
            arcLen = arclength(ctrlX, ctrlY, ctrlZ,'linear');
            samplingFactor = round(arcLen/minD);
            [pt,diffVector] = interparc(samplingFactor, ctrlX, ctrlY, ctrlZ,'pchip');
            if par.r > 0
                nOnRing = 2;
                theta = [0 pi];
                axes = zeros(nOnRing, 3, [size(diffVector,1)]);
                distances=pdist2([ctrlX ctrlY ctrlZ],pt,'euclidean');
                [local_minD,local_idx] = min(distances,[],2);
                ptS=zeros(1,size(diffVector,1));
                for i_seg=1:(size(local_idx,1)-1)
                    ptS(1,local_idx(i_seg):local_idx(i_seg+1))=ctrlS(i_seg);
                end
                for k = 1:size(diffVector,1)
                    [x,y] = pol2cart(theta+(ptS(k)),par.r);%pi*/180
                    xyRing = [zeros(size(x')) x' y'];
                    T = getTransmat(diffVector(k,:),pt(k,:));
                    temp = T*[xyRing.'; ones(1,size(xyRing,1))];
                    axes(:,:,k)= temp(1:end-1,:)';
                end
                axes = permute(axes, [1 3 2]);
                axes = reshape(axes, [],3);
                items(2).XData = axes(:,1);
                items(2).YData = axes(:,2);
                items(2).ZData = axes(:,3);
            else
                items(2).XData = pt(:,1);
                items(2).YData = pt(:,2);
                items(2).ZData = pt(:,3);
            end
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
            [ctrlX,ctrlY,ctrlZ,ctrlS] = getCtrlPointsXYZS(par, obj.internalSettings.numOfCtrlPointSet);
            
            if isempty(obj.ParentObject.locsPrecFactor)
                locsPrecFactor = 1;
            else
                locsPrecFactor = obj.ParentObject.locsPrecFactor;
            end
            
            minD = locsPrecFactor*1; %%0.75
            
            arcLen = arclength(ctrlX, ctrlY, ctrlZ,'linear');
            samplingFactor = round(arcLen/minD);
            [pt,dudt,funthd] = interparc(samplingFactor, ctrlX, ctrlY, ctrlZ,'pchip');
            
            [L,R,K] = curvature(pt);
            derivedPars.curvature = 1./R;
            derivedPars.ctrlpoints.ctrlX = ctrlX;
            derivedPars.ctrlpoints.ctrlY = ctrlY;
            derivedPars.ctrlpoints.ctrlZ = ctrlZ;
            derivedPars.ctrlpoints.ctrlS = ctrlS;
            derivedPars.pt = pt;
            derivedPars.dudt = dudt;
            derivedPars.funthd = funthd;
            derivedPars.curvatureVector = K;
            
            derivedPars.avgCurvature = mean(derivedPars.curvature,'omitnan');
            derivedPars.p75Curvature = prctile(derivedPars.curvature, 75);
            derivedPars.axesRot=ctrlS*180/pi;
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
                lx = startsWith(obj.name,'cx');
                if any(lx)
                    ctrlSetID = str2double(regexprep(obj.name(lx),'\D',''));
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
                    %obj.ParentObject.ParentObject.resetInit;
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
                 % IC220907 addition
                if isfield(obj.ParentObject.ParentObject.parsInit, 'init_locked')
                    obj.ParentObject.ParentObject.lockInit;
                end
            end
        end
    end
end

function struct = defaultParsArgNewCtrlPoint(number)
struct.name = {['cx' num2str(number)],['cy' num2str(number)],['cz' num2str(number)],['s' num2str(number)]};
struct.fix = [0 0 0 0] ;
struct.value = [100*(number-1) 0 0 0];
struct.lb = [-50 -50 -50 -pi];
struct.ub = [50 50 50 pi];
struct.min = [-5000 -5000 -5000 -pi];
struct.max = [5000 5000 5000 inf];
end

function T = getRotMat(p1,t)
p0 = [1 0 0];
p1 = p1./sqrt(p1*p1');
C = cross(p0', p1') ;
D = -dot(p0', p1') ;
NP0 = norm(p0) ; % used for scaling
if ~all(C==0) % check for colinearity
    Z = [0 -C(3) C(2); C(3) 0 -C(1); -C(2) C(1) 0] ;
    R = (eye(3) + Z + Z^2 * (1-D)/(norm(C)^2)) / NP0^2 ; % rotation matrix
else
    %     R = sign(D) * (norm(p1) / NP0) ; % orientation and scaling
    R = eye(3);
end
T = [[R;[0 0 0]] [t'; 1]];
end


function T=getTransmat(r,shift)%,idx,mth
axis=[1 0 0];
r=r/norm(r);
angle= -acosd(dot(r,axis));%- means rotate clockwise!! %atan2d(norm(cross(r,axis)),dot(r,axis)) %LINK missing to matlab and wiki page explaining this
u=cross(axis,r);
if norm(u)==0
    u=u;
    % elseif dot(r,axis)==0
    % angle=angle+1;
    rotMat=eye(3);
    if dot(axis,r)<0
        rotMat(1,1)=-1;
    end
else
    u=u/norm(u);
    rotMat=eye(3)*cosd(angle(1))+sind(angle(1))*[0,-u(3), u(2);u(3),0,-u(1);-u(2),u(1) 0]+(1-cosd(angle(1)))*  [u(1)^2, u(1)*u(2), u(1)*u(3); u(1)*u(2), u(2)^2, u(2)*u(3); u(1)*u(3), u(2)*u(3), u(3)^2 ];

end
T=[[rotMat;[0 0 0]] [shift'; 1]];
end

function [ctrlX,ctrlY,ctrlZ,ctrlS] = getCtrlPointsXYZS(par, numOfCtrlPointSet)
% control points:
      
    for k = 1:numOfCtrlPointSet
        ctrlX(k,:) = par.(['cx' num2str(k)]); % ->
        ctrlY(k,:) = par.(['cy' num2str(k)]); % ->
        ctrlZ(k,:) = par.(['cz' num2str(k)]); % ->
        ctrlS(k,:) = par.(['s' num2str(k)]);
    end
    
end