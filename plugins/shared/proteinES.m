classdef proteinES < handle
    % You can create an object to represent a protien at endocytic sites.
    % You can easily create objects of this class using the segment plug-in
    % proteinES_ctrl.
    % Time-related parameters can be defined in timePar.
    % You can specify the endocytic membrane by modifying the property
    % endocyticPM, which is a basicUnit object.

    properties
        proteinID
        proteinName
        endocyticPM
        capShape
        bottomShape
        timePar
        timeData
        velocityEPM
        zmt0
    end
    methods
        function obj=proteinES(proteinID, proteinName, endocyticPM, capShape, bottomShape, velocityEPM, zmt0)
            obj.proteinID = proteinID;
            obj.proteinName = proteinName;
            obj.endocyticPM = endocyticPM;
            obj.capShape = capShape;
            obj.bottomShape = bottomShape;
            timePar = proteinES.emptyTimePar();
            obj.timePar = timePar;
            obj.timeData = [];
            obj.velocityEPM = velocityEPM;
            obj.zmt0 = zmt0;
        end
        function addTimePoint(obj, aTable)
            obj.timePar = [obj.timePar;aTable];
        end
        function [img, nm] = getTimeImg(obj, t, viewType, imageSize)
            if ~isfield(obj.timeData, proteinES.time2label(t))
                obj.timeData.(proteinES.time2label(t)) = obj.getTime(t);
            end    
            switch viewType
                case 'top'
                    % accumulate all the even points in the space to z=1
                    idx = int16([obj.timeData.(proteinES.time2label(t)).x obj.timeData.(proteinES.time2label(t)).y])+imageSize/2;
                    idxZ = obj.timeData.(proteinES.time2label(t)).z;
                    proImg = accumarray(idx, idxZ,[],@(x) length(x));
                    proImg = [proImg zeros(min(size(proImg)),imageSize-min(size(proImg)))];
                    img = [proImg; zeros(imageSize-min(size(proImg)),imageSize)];
                    img = uint8(round(mat2gray(img)*255));
                    figure(999)
                    image(img, 'CDataMapping', 'scaled');
                    % imwrite(uint8(round(mat2gray(proImg)*255)), [p.folderPath '\' p.imgPath])
                case 'side'
                    % accumulate all the even points in the space to z=1
                    idx = int16([obj.timeData.(proteinES.time2label(t)).x+imageSize/2 obj.timeData.(proteinES.time2label(t)).z+1]);
                    idxZ = obj.timeData.(proteinES.time2label(t)).y+imageSize/2;
                    proImg = accumarray(idx, idxZ,[],@(x) length(x));
                    Size = size(proImg);
                    proImg = [proImg; zeros(imageSize-Size(1),Size(2))];
                    Size = size(proImg);
                    img = [proImg zeros(imageSize, imageSize-Size(2))];
                    img = uint8(round(mat2gray(img)*255));
                    figure(999)
                    image(img', 'CDataMapping', 'scaled');
                    % imwrite(uint8(round(mat2gray(proImg')*255)), [p.folderPath '\' p.imgPath])
            end
            if ~ismember(t, obj.timePar.time)
                thisTP = obj.interpTimePar(t);
                nm = thisTP.nm;
            else
                thisTP = obj.timePar(obj.timePar.time==t,:);
                nm = thisTP.nm;
            end
        end
        function time = getTime(obj, t)
            if t<=3 % determining the last time point that a pre-invagination model should be used.
                if t >= 0
                    thisTP = obj.timePar(obj.timePar.time==999,:);
                else
                    if ismember(t, obj.timePar.time)
                        thisTP = obj.timePar(obj.timePar.time==t,:);
                    else
                        thisTP = obj.interpTimePar(t);
                    end
                end
                % basicUnit(capDepth, capShape, bottomDepth, bottomShape, radius)
                pseudoInner = basicUnit(0, obj.endocyticPM.capShape, obj.endocyticPM.bottomDepth, obj.endocyticPM.bottomShape, obj.endocyticPM.radius+thisTP.shi);
                pseudoOuter = basicUnit(0, obj.endocyticPM.capShape, obj.endocyticPM.bottomDepth, obj.endocyticPM.bottomShape, obj.endocyticPM.radius+thisTP.shi+thisTP.th);
                prePx = fullPixel(pseudoOuter);
                inPseudoOuter = basicUnit.pointsInCy(pseudoOuter.basicComponent3D.x, pseudoOuter.basicComponent3D.y, pseudoOuter.basicComponent3D.z, prePx.x(:), prePx.y(:), prePx.z(:));
                prePx = prePx(inPseudoOuter,:);
                if pseudoInner.radius>0
                    inPseudoInner = basicUnit.pointsInCy(pseudoInner.basicComponent3D.x, pseudoInner.basicComponent3D.y, pseudoInner.basicComponent3D.z, prePx.x(:), prePx.y(:), prePx.z(:));
                    prePx = prePx(~inPseudoInner,:);
                end
                prePx = prePx(prePx.z <= thisTP.zm+thisTP.de/2 & prePx.z >= thisTP.zm-thisTP.de/2 ,:);
            else
                prePx = table([],[],[], 'VariableNames',{'x', 'y', 'z'});
            end
            prePx = prePx(prePx.z>=0,:);
            
            if t>0
                if ismember(t, obj.timePar.time)
                    thisTP = obj.timePar(obj.timePar.time==t,:);
                else
                    thisTP = obj.interpTimePar(t);
                end
                inner = basicUnit(obj.endocyticPM.capDepth+thisTP.shi, obj.endocyticPM.capShape, obj.endocyticPM.bottomDepth, obj.endocyticPM.bottomShape, obj.endocyticPM.radius+thisTP.shi);
                if strcmp(obj.endocyticPM.capShape, obj.capShape) && strcmp(obj.endocyticPM.bottomShape, obj.bottomShape) % Check the shapes of inner and outer columns are the same or not
                    outer = basicUnit(obj.endocyticPM.capDepth+thisTP.shi+thisTP.th, obj.endocyticPM.capShape, obj.endocyticPM.bottomDepth, obj.endocyticPM.bottomShape, obj.endocyticPM.radius+thisTP.shi+thisTP.th);
                end
                Px = fullPixel(outer);
                inOuter = basicUnit.pointsInCy(outer.basicComponent3D.x, outer.basicComponent3D.y, outer.basicComponent3D.z, Px.x(:), Px.y(:), Px.z(:));
                Px = Px(inOuter,:);
                inInner = basicUnit.pointsInCy(inner.basicComponent3D.x, inner.basicComponent3D.y, inner.basicComponent3D.z, Px.x(:), Px.y(:), Px.z(:));
                Px = Px(~inInner,:);
                
                Px = Px(Px.z <= thisTP.zm+thisTP.de/2 & Px.z >= thisTP.zm-thisTP.de/2 ,:);
                
                Px.z = Px.z+ (obj.zmt0 + obj.velocityEPM*t);
            
                innerESC = inner.convert2EScoordinates(obj.zmt0 + obj.velocityEPM*t);
                inInner2 = basicUnit.pointsInCy(innerESC.x, innerESC.y, innerESC.z, prePx.x(:), prePx.y(:), prePx.z(:));
                prePx = prePx(~inInner2,:);
            else
                Px = table([],[],[], 'VariableNames',{'x', 'y', 'z'});
            end
            Px = Px(Px.z>=0,:);
            time = [prePx;Px];
            time = unique(time);
        end
        
        function theAns = getKeyRef(obj, qth, qshi, whichOne, varargin)
            inner = basicUnit(obj.endocyticPM.capDepth+qshi, obj.endocyticPM.capShape, obj.endocyticPM.bottomDepth, obj.endocyticPM.bottomShape, obj.endocyticPM.radius+qshi);
            outer = basicUnit(obj.endocyticPM.capDepth+qshi+qth, obj.endocyticPM.capShape, obj.endocyticPM.bottomDepth, obj.endocyticPM.bottomShape, obj.endocyticPM.radius+qshi+qth);
            switch whichOne
                case 'tip'
                    theAns = outer.totalDepth;
                case 'innerTip'
                    theAns = inner.totalDepth;
                case 'capOri'
                    theAns = inner.bottomDepth;
                case 'outerRad'
                    theAns = outer.radius;
                case 'innerRad'
                    theAns = inner.radius;
                case 'outerRadSpe'
                    idx = outer.basicComponent.y == varargin{1};
                    theAns = outer.basicComponent.x(idx);
            end
        end
        function tp = interpTimePar(obj, t)
            tp = proteinES.emptyTimePar();
            for i=2:6
                thisVal = interp1(obj.timePar.time, obj.timePar{:,i}, t);
                tp{1,i}=thisVal;
            end
            tp{1,1}=t;
            tp.sha={'bar'};
        end
    end
    methods(Static)
        function label = time2label(t)
            label = regexprep(string(t), '\.','d');
            label = regexprep(string(label), '\-','m');
            label = char(strcat('t', label));
        end
        function timePar = emptyTimePar()
            timePar = table([],[],[],[],[],[],[],'VariableNames',{'time', 'th', 'de', 'zm', 'shi', 'nm', 'sha'});
        end

    end
end


