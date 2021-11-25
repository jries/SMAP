classdef NPCPointModel_flexible
    properties
        name = {'ringDistance', 'azimuthalShift', 'radius'};
        fix = [0 0 1] ;
        value = [0 0 53.7];
        lb = [0 -inf 0];
        ub = [70 inf 0];
        min = [0 -inf 0];
        max = [70 inf 100];
        modelType = 'discrete'
        dimension = 3;
    end
    methods
        function obj = NPCPointModel_flexible
        end
    end
    methods (Static)
        function [model, p]= reference(par, dx)
        % set additional parameters of the model
        copyPerCorner = 2;
        cornerNum = 8;
        useSecondRing = true;
        cornerRange = 12;
        
        % corner:
        cornerPos = linspace(0,2*pi,cornerNum+1)-pi/cornerNum;
        cornerPos = cornerPos(1:end-1);
        
        % copies in one corner:
        copyPos = linspace(0,cornerRange*pi/180,copyPerCorner);
        
        % copies in all corners of one ring
        allCopiesTheta = copyPos + cornerPos';
        
        % vecterization
        fn = fieldnames(par);
        maxParElem = 0;
        for k = 1:length(fn)
            maxParElem = max([maxParElem length(par.(fn{k}))]);
        end
        tempTheta = zeros([1 maxParElem]);              % this is for the vecterization
        allCopiesTheta  = allCopiesTheta(:)+tempTheta;
        par.radius = tempTheta+par.radius;
        
        % assign the radius for each copy
        allCopiesRho = repelem(par.radius', size(allCopiesTheta,1));
        
        % convert from polar coordinates
        [allCopiesX,allCopiesY]= pol2cart(allCopiesTheta(:), allCopiesRho(:));
        allCopiesX = reshape(allCopiesX,size(allCopiesTheta));
        allCopiesY = reshape(allCopiesY,size(allCopiesTheta));
        
        % for the second ring
        if useSecondRing
            allCopiesThetaRing2 = allCopiesTheta+par.azimuthalShift.*pi./180;
            [allCopiesXRing2,allCopiesYRing2]= pol2cart(allCopiesThetaRing2(:), allCopiesRho(:));
            allCopiesXRing2 = reshape(allCopiesXRing2,size(allCopiesTheta));
            allCopiesYRing2 = reshape(allCopiesYRing2,size(allCopiesTheta));
        end
        halfRingDist = par.ringDistance./2;
        
        allCopiesZ = zeros(size(allCopiesX))+halfRingDist;
        allCopiesZRing2 = zeros(size(allCopiesXRing2))-halfRingDist;
        model.x = [allCopiesX;allCopiesXRing2];
        model.y = [allCopiesY;allCopiesYRing2];
        model.z = [allCopiesZ;allCopiesZRing2];
        model.channel = ones(size(model.x));
        model.n = ones(size(model.x));

        %     likelihood = gaussDist(model,double(locs.ynmrot), double(locs.xnmrot), double(locs.znm), 8,8,10);
        %     likelihoodBoun = gaussDist(model,0,0,0, 16,16,20);

        p.cornerRange=cornerRange;
        p.copyPerCorner=copyPerCorner;
        p.cornerNum=cornerNum;
        end
    end
end
