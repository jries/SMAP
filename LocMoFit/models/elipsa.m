classdef elipsa<geometricModel
    
    % log
    %   - 201229: change the sign of the ring twist
    methods
        function obj = elipsa(varargin)
            obj@geometricModel(varargin{:});
            % Define parameters that can be altered during fitting here:
            obj.name = {'a', 'b'}; % parameter names
            obj.fix = [0 0] ;                                                    % fix to a constant or not
            obj.value = [45 65];                                               % initial guess
            obj.lb = [-5 -5];                                                   % relative lower bound
            obj.ub = [5 5];                                                % relative upper bound
            obj.min = [35 55];                                               % absolute lower bound
            obj.max = [55 75];                                              % absolute upper bound
            
            % Define discrite parameters here:
            %obj.internalSettings.copyPerCorner = 2;
            obj.internalSettings.cornerNum = 24;
            %obj.internalSettings.useSecondRing = true;
            
            % Define other properties here:
            obj.modelType = 'discrete';
            obj.modelTypeOption = {'discrete'};
            obj.dimension = 3;
            
        end
        
        function [model, p]= reference(obj, par, dx)
        % Sample coordinates of the model as reference.
        % --- Syntax ---
        % [model, p]= reference(obj, par, dx)
        % --- Arguments ---
        % -- Input --
        % obj:
        % par: a structure object. Its fieldnames should be the names of
        % parameters, and their correspoinding content should be the
        % parameter values.
        % dx: sampling rate.
        % -- Output --
        % model: a structure object. Its fieldnames should be x, y, z, and
        % n, indicating the xyz position amplitude n of the sampled model
        % points.
        % p: additional information of the model.
        
        
        % set additional parameters of the model
        ip = inputParser;
        fn = fieldnames(obj.internalSettings);
        for k = 1:length(fn)
            ip.addParameter(fn{k}, obj.internalSettings.(fn{k}));
        end
        a = {};
        parse(ip,a{:})
        pResults = ip.Results;
        
        %copyPerCorner = pResults.copyPerCorner;
        cornerNum = pResults.cornerNum;
        %useSecondRing = pResults.useSecondRing;
        %cornerRange = par.cornerDegree;
        
        % corner:
%         cornerPos = linspace(0,2*pi,cornerNum+1);
%         cornerPos = cornerPos(1:end-1);
        
%         % copies in one corner:
%         for k=length(cornerRange):-1:1
%             copyPos(:,k) = linspace(0,cornerRange(k)*pi/180,copyPerCorner)-cornerRange(k)*pi/360; % the last term is for aligning the corners among different cornerRange
%         end
%         
%         % copies in all corners of one ring
%         cornerPosNew(1,1,:) = cornerPos;
%         allCopiesTheta = copyPos + cornerPosNew;
        
        if isempty(obj.ParentObject.locsPrecFactor)
            locsPrecFactor = 1;
        else
            locsPrecFactor = obj.ParentObject.locsPrecFactor;
        end
        %minDist = locsPrecFactor*dx;      % change this if you need more points

        % corner:
        theta = linspace(0,2*pi,cornerNum+1);
        theta = theta(1:end-1);
        
        
        oneRingX = par.a*cos(theta);     
        oneRingY= par.b*sin(theta)  ;
        

        model.x = oneRingX.';
        model.y = oneRingY.';
        model.z = zeros(length(oneRingX),1);

        model.channel = ones(size(model.x));
        model.n = ones(size(model.x));
        p = [];

        %     likelihood = gaussDist(model,double(locs.ynmrot), double(locs.xnmrot), double(locs.znm), 8,8,10);
        %     likelihoodBoun = gaussDist(model,0,0,0, 16,16,20);

        
        %p.cornerRange=cornerRange;
        %p.copyPerCorner=copyPerCorner;
        p.cornerNum=cornerNum;
        end
        function derivedPars = getDerivedPars(obj, pars)
            derivedPars = [];
        end
    end
end
