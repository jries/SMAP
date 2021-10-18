classdef LocMoFit<matlab.mixin.Copyable
    % This is the class for fitting a set of geometric models to SMLM data.
    %
    % This class help user to set up and handle the model fitting.
    %
    % If you would like to perform multi-step fitting, please create one
    % SMLMModelFit object for each step.
    properties
        dimension               % ??? The dimension of the data.
        allParsArg              % All arguments of the parameters.
        parsInit                % ???
        fitInfo                 % Additional information acquired by fitting.
        externalInfo = []       % External info, which were not acquired by the fit but rather used for alignment or so.
        model                   % The model objects.
        solver                  % The solver options.
        converterRules          % The defined conversion rules.
        converterUserDefined    % The temporary variables that the user defined.
        roiSize = 300;          % The size of region of interest in unit of nm.
        imgExtension = 0;       % The extention of image size than the roi size.
        dataDim = 3;            % The dimension of the data.
        allModelLayer           % The unique numbers of layers used by models.
        modelLayer              % The layer for the corresponding model.
        sigmaCascade = 1;       % The sigma factor for each step of the sequential cascade fit.
        modelVerCascade = 1;    % For locs model only. Allow the user to specify different versions of the same model for different cascade steps.
        refPoint_spacing = 0.75;% The spacing bwteen sampled ref points. In the unit of sigma.
        numOfLocsPerLayer       % The numbers of localizations per layer.
        representiveLocprec     % The representative localization precision.
        display                 % Settings for displaying figures.
        lutLayer
        objFunType = 'likelihood'; % The type of the objective function.
        useCompensation = true;
        compensationFactor
        weightLayer = 1;       % The weights of layers.
        fitterInfo
        status = 'finished';   % can be 'initial', 'iterative' or 'finished'
        advanceSetting          % Advanced settings are defined here.
    end
    properties (Transient)
        linkedGUI               % If GUI is used, it will be saved here.
        converterSource         % The information source for the converter.
        locs
        currentFitPars         % the current parameters of the model being evaluated.
        handles                % handles to graphic objects
        temp                    % all temp info is saved here.
    end
    properties (Dependent)
        numOfModel              % The number of models used in this fitting step.
        numOfLayer              % The number of layers.
        roiAreaVol              % The area/volume of the roi;
    end
    properties (Hidden)
        supportedSolver = {'fmincon', 'fminsearchbnd','particleswarm'};
        currentCascadeStep = 1; % Indicating the curent cascasde step.
    end
    methods
        function obj = LocMoFit(varargin)
            % Construct the object of the class 'LocMoFit'
            %
            % Usage:
            %   obj = LocMoFit(Name-value)
            %
            % Args:
            %   Name-value pairs:
            %       * 'SolverName': one of 'fmincon', 'fminsearchbnd', and 'particleswarm'.
            %       * 'SolverOptions': 
            %       * 'DataDim': 
            %       * 'TestLocs': 
            % Returns:
            %   obj: an LocMoFit object.
            %
            % NOTE:
            %   Please create an :class:`LocMoFit` object for each step of
            %   fitting.
            %
            % See also:
            %   :func:`setModel`
            %
            
            locs = [];
            p = inputParser;
            p.addParameter('SolverName', 'fminsearchbnd')
            p.addParameter('SolverOptions', {})
            p.addParameter('DataDim', 2)
            p.addParameter('TestLocs',locs);
            p.parse(varargin{:})
            
            pArgin = p.Results;
            obj.solver.SolverName = pArgin.SolverName;
            obj.solver.SolverOptions = pArgin.SolverOptions;
            obj.fitterInfo.possibleSolver = {'fminsearchbnd','particleswarm','fmincon'};
            if pArgin.DataDim==3||isfield(pArgin.TestLocs, 'znm')
                obj.dataDim = 3;
            else
                obj.dataDim = 2;
            end

            obj.loadListener;
            obj.initAdvanceSetting;
        end
        function loadListener(obj)
            addlistener(obj, 'mParsArgModified', @mParsArgModified_callback);
        end
        
        %% separate functions
        converter(obj, source, rule, target)
        [totalIntensity,wk,LLctrl] = intensityCal(obj, fitPars, locs)
        [ax,finalImg] = plot(obj,locs,varargin)
        modPoint = getModPoint(obj, modelSamplingFactor)
        layerPoint = getLayerPoint(obj, modSamplingF)
        modPoint = transModPoint(obj, modPoint)
        fit(obj, locs, varargin)
        convertNow(obj, locs)
        stop = plotFreeRot(obj, varargin)
        fig = plotFixRot(obj, varargin)
        locs = locsHandler(obj, locs, lParsVal,modelID,varargin)
        ax = rotCoordNMkImg(obj, varargin)
        [lPars, mPars,fitPars] = inputPar2struct(obj, k, fitPars, lPars, mPars, offset)
        stop = optimoPlotMod(obj,x,optimValues,state,varargin)
        item = getThings2Plot(obj,varargin);
        stop = optimOutputPar(obj,x,optimValues,state,varargin)
        showFitResult(obj,varargin);
        
        d = distance2RefPoint(obj, refPoint,varargin)
        %% model related functions
        function setModel(obj,model,modelId)
            % Adding one single model to the LocMoFit object according to the modelId.
            % Initiation of all arguments of the parameters (allParsArg).
            %
            % Usage:
            %   setModel(obj,model,modelId)
            %
            % Args:
            %   model: an SMLMModel object or sub-object.
            %   modelId: the ID of the model being added.
            
            % add parent to the model obj and mount the model
            model.addParent(obj);
            obj.model{modelId} = model;
            obj.model{modelId}.ID = modelId;
            lPars = defaultLpars(obj, model.dimension);
            if model.dimension~=obj.dataDim
                switch obj.dataDim
                    case 2
                        error('Invalid: You are adding a 3D model for fitting 2D data.')
                        return
                    case 3
                        obj.dataDim = 2;
                        warning('You are adding a 2D model for fitting 3D data. Switching to 2D fit.')
                end
            end
            
            fn = fieldnames(lPars);
            
            % check the layer is new or not
            if isempty(obj.allModelLayer)           % if there is no existing layer, that means this is the first model
                obj.allParsArg.model = [];
                obj.allParsArg.type = [];
                obj.allParsArg.label = [];
                for l = 1:length(fn)
                    obj.allParsArg.(fn{l}) = [];
                end
            end
            
            % fetch defaults
                %--for lPars
            obj.allParsArg.model = [obj.allParsArg.model; repelem(modelId,length(lPars.name),1)];
            obj.allParsArg.type = [obj.allParsArg.type; repelem({'lPar'},length(lPars.name),1)];
            obj.allParsArg.label = [obj.allParsArg.label; repelem(cellstr(''),length(lPars.name),1)];
            for l = 1:length(fn)
                obj.allParsArg.(fn{l}) = [obj.allParsArg.(fn{l}); lPars.(fn{l})];
            end
            
                %--for mPars
            mPars = obj.model{modelId}.mPars;
            if ~isempty(mPars)
                for l = 1:length(fn)
                    obj.allParsArg.(fn{l}) = [obj.allParsArg.(fn{l}); mPars.(fn{l})'];
                end
                obj.allParsArg.model = [obj.allParsArg.model; repelem(modelId,length(mPars.name),1)];
                obj.allParsArg.type = [obj.allParsArg.type; repelem({'mPar'},length(mPars.name),1)];
                obj.allParsArg.label = [obj.allParsArg.label; repelem(cellstr(''),length(mPars.name),1)];
            end
            obj.updateLayer;

        end
        
        function addModel(obj,model)
            % Add a specified model to the end of the list of models.
            %
            % Usage:
            %   addModel(obj,model)
            %
            % Args:
            %   obj: an LocMoFit object.
            %   model: an SMLMModel object or sub-object.
            % TODO:
            %   flag1: obj.initLParSelector for will be available when
            %   implemnted
            
            lastMod = length(obj.model);
            if 0
                % TODO: flag1
                obj.initLParSelector(lastMod+1, false);
            end
            obj.setModel(model,lastMod+1);
            obj.model{lastMod+1}.ID = lastMod+1;
            obj.updateAProp('numOfLayer');
        end
        
        function changeModel(obj, newModel, modelNumber)
            % Remove the old corresponding parameters and add a new model to overwrite the old model with the same ID.
           	% 
            % Usage:
            %   changeModel(obj, newModel, modelNumber)
            %
            % Args:
            %   obj: an LocMoFit object.
            %   newModel: an SMLMModel object or sub-object being added.
            %   modelNumberRemove: the ID of the model being added.
            
            lModel = ismember(obj.allParsArg.model, modelNumber);
            fn = fieldnames(obj.allParsArg);
            
            % Set parameters of the model to null
            for k = 1:length(fn)
                obj.allParsArg.(fn{k})(lModel) = [];
            end
            obj.setModel(newModel,modelNumber);
        end
        
        function updateModel(obj, modelNumber)
            % Respond to the change of any internal setting of the geometric model.
            
            % remove the corresponding parameters
            lModel = ismember(obj.allParsArg.model, modelNumber);
            lMPar = ismember(obj.allParsArg.type, 'mPar');
            fn = fieldnames(obj.allParsArg);
            parsNameNew = obj.model{modelNumber}.mPars.name;
            
            % Detect the parameters being removed
            lKept = zeros(size(obj.allParsArg.name));
            for k = 1:length(parsNameNew)
                lKept(lModel&lMPar) = lKept(lModel&lMPar) + strcmp(obj.allParsArg.name(lModel&lMPar), parsNameNew{k});
            end
            for k = 1:length(fn)
                obj.allParsArg.(fn{k})(lModel&lMPar&~lKept) = [];
            end
            
            % Detect the parameters being added
            lModel = ismember(obj.allParsArg.model, modelNumber);
            lMPar = ismember(obj.allParsArg.type, 'mPar');
            parsNameOld = obj.allParsArg.name(lModel&lMPar);
            lOld = zeros(size(parsNameNew));
            for k = 1:length(parsNameOld)
                lOld = lOld + strcmp(parsNameNew, parsNameOld{k});
            end
            indNew = find(~lOld);
            fn = {'fix','value','lb','ub','min','max'};
            for l = 1:length(indNew)
                for k = 1:length(fn)
                    obj.allParsArg.(fn{k})(end+1) = obj.model{modelNumber}.mPars.(fn{k})(indNew(l));
                end
                obj.allParsArg.name{end+1} = parsNameNew{indNew(l)};
                obj.allParsArg.model(end+1) = modelNumber;
                obj.allParsArg.label{end+1} = char();
                obj.allParsArg.type{end+1} = 'mPar';
            end
            
            lOffSet = strcmp(obj.allParsArg.type, 'offset');
            indOffSet = find(lOffSet);
            for k = 1:length(fn)
                obj.allParsArg.(fn{k})(end+1) = obj.allParsArg.(fn{k})(lOffSet);
                obj.allParsArg.(fn{k})(indOffSet) = [];
            end
            obj.allParsArg.name{end+1} = obj.allParsArg.name{lOffSet};
            obj.allParsArg.name(indOffSet) = [];
            obj.allParsArg.model(end+1) = obj.allParsArg.model(lOffSet);
            obj.allParsArg.model(indOffSet) = [];
            obj.allParsArg.label{end+1} = obj.allParsArg.label{lOffSet};
            obj.allParsArg.label(indOffSet) = [];
            obj.allParsArg.type{end+1} = obj.allParsArg.type{lOffSet};
            obj.allParsArg.type(indOffSet) = [];
        end
        
        function img = getImage(obj,modelNumber,varargin)
            % Get an image of the specified model at the best parameters.
            mPars = obj.exportPars(modelNumber, 'mPar');
            img = obj.model{modelNumber}.getImage(mPars,varargin{:});
        end
        
        %% for layers
        function updateLayer(obj)
            % Manage layer-dependent offsets.
            obj.modelLayer = getFieldAsVector(obj.model,'layer'); % the layers where models belong to
            oldLayerInUse = obj.allModelLayer;
            obj.allModelLayer = unique(obj.modelLayer);
            
            % Remove layers not in use.
            for k = 1:length(oldLayerInUse)
                lLayerRm = ~ismember(oldLayerInUse(k),obj.allModelLayer);
                if lLayerRm
                    obj.rmPar(['m9' num2str(oldLayerInUse(k)) '.offset.weight'])
                end
            end
            
            % Add new layers.
            indNewLayer = find(~ismember(obj.allModelLayer, oldLayerInUse));
            if ~isempty(indNewLayer)
                offset = defaultOffset(obj, {'weight'});
                fnO = fieldnames(offset);
                for k = 1:length(indNewLayer)
                    modelNumber = 9*10 + obj.allModelLayer(indNewLayer(k));
                    obj.initLParSelector(modelNumber, true);
                    obj.allParsArg.model = [obj.allParsArg.model; repelem(modelNumber,length(offset.name),1)];
                    obj.allParsArg.type = [obj.allParsArg.type; repelem({'offset'},length(offset.name),1)];
                    obj.allParsArg.label = [obj.allParsArg.label; repelem(cellstr(''),length(offset.name),1)];
                    for l = 1:length(fnO)
                        obj.allParsArg.(fnO{l}) = [obj.allParsArg.(fnO{l}); offset.(fnO{l})];
                    end
                end
                
            end
            
            % Reorder the pars to make sure offsets' always stay at the end
            lReorder = find(strcmp(obj.allParsArg.type, 'offset'));
            fnReorder = fieldnames(obj.allParsArg);
            for k = 1:length(fnReorder)
                obj.allParsArg.(fnReorder{k})(end+1:end+length(lReorder)) = obj.allParsArg.(fnReorder{k})(lReorder);
                obj.allParsArg.(fnReorder{k})(lReorder) = [];
            end
            obj.allParsArg.fix = logical(obj.allParsArg.fix);
        end
        
        function updateAProp(obj, oneProp)
            switch oneProp
                case 'numOfLayer'
                    
                otherwise
                    warning('The property is not supported.')
            end 
        end
        
%         function flagNewLayer = changeLayer(obj, newLayer)
%             % Check whether the change of layer introduce a new layer.
%             %
%             % --- Syntax ---
%             % flagNewLayer = changeLayer(obj, newLayer)
%             % --- Discription ---
%             % flagNewLayer: logical. Indicating the layer is new or not.
%             % newLayer: an integer indicating the new layer number.
%             
%             % First, always assue the layer has been added.
%             flagNewLayer = 0;
%             
%             % Second, compare the new layer to the existing layers. If the
%             % new layer is really new, then change the flag.
%             if ~ismember(newLayer, obj.allModelLayer)  
%                 obj.allModelLayer(end+1) = obj.model{end}.layer;
%                 obj.numOfLayer = obj.numOfLayer+1;
%                 flagNewLayer = 1;
%             end
%         end
        
        %% Parameter related methods
        function parsArg = subParsArg(obj, model, varargin)
            p = inputParser;
            p.addParameter('Type', []);
            parse(p, varargin{:})
            results = p.Results;
           
            if isempty(model)% for layer settings
                lModel = obj.allParsArg.model > 90;
            else
                lModel = obj.allParsArg.model == model;
            end
            
            if isempty(results.Type)
                lType = true(size(obj.allParsArg.type));
            else
                lType = strcmp(obj.allParsArg.type, type);
            end
            fn = fieldnames(obj.allParsArg);
            for k = 1:length(fn)
                parsArg.(fn{k}) = obj.allParsArg.(fn{k})(lModel&lType);
            end
        end
        
        function defineRotPar(obj)
            
        end
        
        function [lb,ub,value, min, max] = prepFit(obj)
            % Reshape arguments of fit pars: lb, ub, init(value), and min/max
            
            % Simply subsetting the arguments for fit parameters.
            indFit = ~obj.allParsArg.fix;
            % Convert the arguments into multiple vectors.
            lb = obj.allParsArg.lb(indFit);
            ub = obj.allParsArg.ub(indFit);
            value = obj.allParsArg.value(indFit);
            min = obj.allParsArg.min(indFit);
            max = obj.allParsArg.max(indFit);
        end
        
        function k = numOfFreePar(obj)
            indFit = ~obj.allParsArg.fix;
            k = sum(indFit);
        end
        
        function updatePars(obj,newParsArg)
            fn = fieldnames(newParsArg);
            for k=1:length(newParsArg.name)
                indNames = ismember(obj.allParsArg.name, newParsArg.name{k})&...
                    ismember(obj.allParsArg.type, newParsArg.type(k))&...
                    ismember(obj.allParsArg.model, newParsArg.model(k));
                for l = 1:length(fn)
                    obj.allParsArg.(fn{l})(indNames) = newParsArg.(fn{l})(k);
                end
            end
            indFix = obj.allParsArg.fix;
            obj.allParsArg.value(indFix) = obj.allParsArg.value(indFix);
        end
        function allParsArgTable = viewPars(obj)
            % Display the parsArg as a table.
            localAllParsArg = obj.allParsArg;
            
            % Don't display the initial guess.
            if isfield(localAllParsArg, 'init')
                localAllParsArg = rmfield(localAllParsArg, 'init');
            end
            allParsArgTable = struct2table(localAllParsArg);
        end
        function pars = exportPars(obj, modelID, type)
            % Get parameters for a specific model.
            %
            % --- Syntax ---
            % pars = exportPars(obj, modelID, type)
            % 
            % --- Description ---
            % pars: a structral array with parameter names as field names.
            % obj: an LocMoFit object.
            % modelID: the ID of the model where you want to get parameters.
            % type: either 'lPar', 'mPar' or 'allPar', specifying the type of parameters you want to get.
            
            
            indModel = obj.allParsArg.model==modelID;
            if strcmp(type, 'allPar')
                indType = true(size(obj.allParsArg.type));
            else
                indType = ismember(obj.allParsArg.type,type);
            end
            fn = obj.allParsArg.name(indModel&indType);
            values = obj.allParsArg.value(indModel&indType);
            for k = 1:length(fn)
                pars.(fn{k}) = values(k);
            end
            if length(fn) == 0
                pars = [];
            end
        end
        
        function [val,ind] = wherePar(obj, parId)
            % [Replaced] see getVariable().
            % 200731: this function has been replaced getVariable()
            disp('[Obselete] wherePar() will be replaced by getVariable().')
            [val,ind] = getVariable(obj, parId);
        end
        
        function [val,ind] = getVariable(obj, ID)
            % Search for a variables in where info is potentially stored 
            % and report its location and value.
            % ID shold look like par.m1.lPar.x or directly the variable
            % name.
            parts = strsplit(ID,'.');
            if length(parts) == 4
                model = replace(parts{2},'m','');
                type = parts{3};
                name = parts{4};
                lModel = ismember(obj.allParsArg.model, str2num(model));
                lType = ismember(obj.allParsArg.type, type);
                lName = ismember(obj.allParsArg.name, name);
                ind = find(lModel&lType&lName);
                if ~isempty(ind)
                    val = obj.allParsArg.value(ind);
%                     disp('[Obselete] 4-element IDs might not be supported in the furture version. Please use 2-element IDs instead.')
                else
                    val = [];
                    ind = [];
                end
            elseif length(parts) == 3
                model = replace(parts{1},'m','');
                type = parts{2};
                name = parts{3};
                lModel = ismember(obj.allParsArg.model, str2num(model));
                lType = ismember(obj.allParsArg.type, type);
                lName = ismember(obj.allParsArg.name, name);
                ind = find(lModel&lType&lName);
                if ~isnan(obj.allParsArg.value(ind))
                    val = obj.allParsArg.value(ind);
%                     disp('[Obselete] 4-element IDs might not be supported in the furture version. Please use 2-element IDs instead.')
                else
                    val = [];
                    ind = [];
                end
            elseif length(parts)==2
                if ~isempty(regexp(parts{1}, '^m\d*','once'))
                    model = replace(parts{1},'m','');
%                     ('Check here if there is an error related to background.')
%                     if str2num(model)>90
%                         model = num2str(str2num(model)-90);
%                     end
                    name = parts{2};

                    if sum(ismember(obj.allParsArg.name, name))
                        lModel = ismember(obj.allParsArg.model, str2num(model));
                        lName = ismember(obj.allParsArg.name, name);
                        ind = find(lModel&lName);
                        val = obj.allParsArg.value(ind);
                    elseif isfield(obj.fitInfo.modelPar_internal{str2num(model)}, name)
                        val = obj.fitInfo.modelPar_internal{str2num(model)}.(name);
                        ind = ['.fitInfo.modelPar_internal{' model '}.' name];
                    elseif isfield(obj.fitInfo.derivedPars{str2num(model)}, name)
                        val = obj.fitInfo.derivedPars{str2num(model)}.(name);
                        ind = ['.fitInfo.derivedPars{' model '}.' name];
                    else
                        val = [];
                        ind = [];
                    end
                elseif ~isempty(regexp(parts{2}, '^l\d*','once'))
                    layer = replace(parts{2},'l','');
                    name = parts{1};
                    val = obj.fitInfo.(name)(str2num(layer));
                    ind = ['.fitInfo.' name '(' layer ')'];
                else
                    val = [];
                    ind = [];
                end
            elseif length(parts)==1
                parts = parts{1};
                if isfield(obj.converterUserDefined, parts)
                    val = obj.converterUserDefined.(parts);
                    ind = ['.converterUserDefined.' parts];
                elseif isfield(obj.fitInfo, parts)
                    val = obj.fitInfo.(parts);
                    ind = ['.fitInfo.' parts];
               elseif isfield(obj.fitInfo, 'derivedPars')&&isfield(obj.fitInfo.derivedPars, parts)
                    val = obj.fitInfo.derivedPars.(parts);
                    ind = ['.fitInfo.derivedPars.' parts];
                elseif isfield(obj.externalInfo, parts)
                    val = obj.externalInfo.(parts);
                    ind = ['.externalInfo.' parts];
                else
                    val = [];
                    ind = [];
                end
            else
                val = [];
                ind = [];
            end
            if isempty(val)
                val = nan;
            end
        end
        
        function offset = exportOffset(obj,fitPars)
            offset = {};
            for ch = 1:obj.numOfLayer
                layer = obj.allModelLayer(ch);
                % fetch the weight of the offset
                indLayer = obj.allParsArg.model == 90+layer;
                indOs = ismember(obj.allParsArg.type,'offset');
                indFit = ~obj.allParsArg.fix;
                fn = obj.allParsArg.name(indLayer&indOs&indFit);
                if ~isempty(fn)
                    offset = LocMoFit.pars2struct(fn, offset, fitPars, 90+layer, 'AddUp', 0);
                    fitPars = fitPars(:,length(fn)+1:end);
                else
                    fn = obj.allParsArg.name(indLayer&indOs&~indFit);
                    value = obj.allParsArg.value(indLayer&indOs&~indFit);
                    offset = LocMoFit.pars2struct(fn, offset, value', 90+layer, 'AddUp', 0);
                end
            end
        end
        
        function setParArgBatch(obj, parArg, varargin)
            nPar = length(parArg.name);
            for k = 1:nPar
                parIdFull = ['pars.m' num2str(parArg.model(k)) '.' parArg.type{k} '.' parArg.name{k}];
                parId = ['m' num2str(parArg.model(k)) '.' parArg.type{k} '.' parArg.name{k}];
                if ~isempty(obj.getVariable(parIdFull))
                    obj.setParArg(parId,'value', parArg.value(k), 'lb', parArg.lb(k),'ub', parArg.ub(k),'fix', parArg.fix(k),'min', parArg.min(k),'max', parArg.max(k),'label', parArg.label(k))
                end
            end
        end
        
        
        function setParArg(obj, parId, varargin)
            % parId looks like this: 'm1.lPar.x', where m1 means model
            % 1, and x can be any name of the parameters
            p = inputParser;
            p.addParameter('lb', []);
            p.addParameter('ub', []);
            p.addParameter('value', []);
            p.addParameter('fix', []);
            p.addParameter('label', []);
            p.addParameter('min', []);
            p.addParameter('max', []);
            p.addParameter('Remove', false);
            fn = {'lb','ub','value','fix','label','min','max'};
            parse(p,varargin{:});
            results = p.Results;
            lRm = results.Remove;
            results = rmfield(results, 'Remove');
            % Identify the parameter to modify
            [~,ind] = obj.getVariable(['pars.' parId]);
            
            % Check through the modifiable fields
            for k = 1:length(fn)
                if ~isempty(ind)
                    if ~isempty(results.(fn{k}))||lRm
                        if isequal(fn{k}, 'label')&&~isempty(obj.allParsArg.(fn{k}){ind})
                            obj.allParsArg.(fn{k}){ind} = ['__' results.(fn{k})]; % put 2 underscore before the label to identify that it is a label
                        else
                            obj.allParsArg.(fn{k})(ind) = results.(fn{k});
                        end
                    end
                end
            end
        end
        function addPar(obj, parArg)
            % Add a parameter.
            % Usage:
            %   obj.addPar(parArg)
            % Arg:
            %   *parArg: a cell in the order of 'model','type','name',
            %   'lb','ub','value','fix','label','min' and 'max'
            fn = {'model','type','name', 'lb','ub','value','fix','label','min','max'};
            if iscell(parArg)
                % Check through the modifiable fields
                for k = 1:length(fn)
                    obj.allParsArg.(fn{k})(end+1) = parArg{k};
                end
            else
                % Check through the modifiable fields
                nPar = length(parArg.name);
                for k = 1:length(fn)
                    obj.allParsArg.(fn{k})(end+1:end+nPar) = parArg.(fn{k});
                end
            end
        end
        function rmPar(obj, parId)
            % Remove a parameter
            [~,ind] = obj.getVariable(['pars.' parId]);
            fn = {'model','type','name', 'lb','ub','value','fix','label','min','max'};
            % Check through the modifiable fields
            for k = 1:length(fn)
                if ~isempty(ind)
                   obj.allParsArg.(fn{k})(ind) = [];
                end
            end
        end
        function resetInit(obj)                                   % reset the init using the original values (before fitting)
            if isfield(obj.parsInit, 'init_locked')
                obj.allParsArg.value = obj.parsInit.init_locked;
            elseif isfield(obj.parsInit, 'init')
                obj.allParsArg.value = obj.parsInit.init;
            end
        end
        function saveInit(obj)
            obj.parsInit.init = obj.allParsArg.value;
        end
        
        function [parId,subParsArgTemp] = getAllParId(obj, modelnumber, varargin)
            % Export Id for all parameters given a model number
            % Show all parId when modelnumber is not specified.
            %
            % Usage:
            %   [parId,subParsArgTemp] = obj.getAllParId(modelnumber, varargin)
            %
            % Args:
            %   modelnumber: an LocMoFit object.
            %   Name-value pairs:
            %       'form': either 'short', 'long', 'auxiliary'
            % Returns:
            %   modCoord: reference coordinates.
            
            p = inputParser;
            p.addParameter('type','main');
            p.addParameter('form','short');
            p.parse(varargin{:})
            results = p.Results;
            
            switch results.type
                case 'main'
                    flag = [1 0];
                case 'all'
                    flag = [1 1];
                case 'auxiliary'
                    flag = [0 1];
            end
            
            parId = [];
            if flag(1)
                % Main parameters
                if exist('modelnumber','var')&&~isempty(modelnumber)
                    modelnumberStr = num2str(modelnumber);
                    subParsArg = obj.subParsArg(modelnumber);
                    fn = fieldnames(obj.allParsArg);
                    for k = 1:length(fn)
                        if iscell(subParsArg.(fn{k}))
                            subParsArgTemp.(fn{k}) = char(subParsArg.(fn{k}));
                        else
                            subParsArgTemp.(fn{k}) = subParsArg.(fn{k});
                        end
                    end

                    for l = length(subParsArg.model):-1:1
                        parId{l} = ['m' modelnumberStr '.' subParsArg.type{l} '.' subParsArg.name{l}];
                    end
                else
                    parId = strcat('m', cellstr(string(obj.allParsArg.model)), '.', obj.allParsArg.type, '.', obj.allParsArg.name);
                    subParsArgTemp = obj.allParsArg;
                end 
            end
            
            if flag(2)
                % Auxiliary parameters or info
                layerDependent = {'numOfLocsPerLayer'; 'numOfLocsPerLayer_BGFree'};
                
                fnLayer = [];
                for l = 1:obj.numOfLayer
                    fn = strcat('fitInfo.', layerDependent, ['.l' num2str(l)]);
                    fnLayer = [fnLayer; fn];
                end
                
                fnInternal = [];
                if ~isempty(obj.fitInfo.modelPar_internal)
                    for m = 1:obj.numOfModel
                        fn = fieldnames(obj.fitInfo.modelPar_internal{m});
                        fn = strcat('fitInfo.modelPar_internal.', ['m' num2str(m) '.'], fn);
                        fnInternal = [fnInternal; fn];
                    end
                end
                
                fnDerived = [];
                if ~isempty(obj.fitInfo.derivedPars)
                    for m = 1:obj.numOfModel
                        fn = fieldnames(obj.fitInfo.derivedPars{m});
                        fn = strcat('fitInfo.derivedPars.', ['m' num2str(m) '.'], fn);
                        fnDerived = [fnDerived; fn];
                    end
                end
                parId = [parId; fnInternal; fnDerived; fnLayer];
            end
            if strcmp(results.form, 'short')
                % turn the IDs into short form
                parId = regexprep(parId,'\.[lm]Par\.', '.');
                parId = regexprep(parId,'\.offset.weight$', '.offset');
                parId = regexprep(parId,'^fitInfo\.(\w*)\.(l\d)$', '$1\.$2');
                parId = regexprep(parId,'^(fitInfo\.\w*\.)', '');
            end
        end

        function interalLPars = convert2InteralLPar(obj, lPars)
            % TODO: define the conversion for lPars in the furture
            % Now this is just a placeholder
        	interalLPars = lPars;
        end
        function internalOffset = convert2InteralOffset(obj, offset)
            % This function converts any valid offset to the interanl
            % offset (weight)
            for l = 1:obj.numOfLayer
                oneOffset = offset{90+l};
                fn = fieldnames(oneOffset);
                for k = length(fn):-1:1
                    if strcmp(fn,'density')
                        internalOffset{90+l}.weight = obj.convertBG('weight', l, oneOffset.density);
                    else
                        internalOffset{90+l} = oneOffset;
                    end
                end
            end
        end
        
        function density = getBGDensity(obj)
            % density here is defined as the projected density.
            % this unit is locs/um^2
            numOfLayer = obj.numOfLayer;
            for l = numOfLayer:-1:1
                density(l) = obj.getBG('density', l);
            end
        end
        
        function BG = getBG(obj, queriedForm, layer)
            optionValue = obj.advanceSetting.(['m9' num2str(layer) '_background']).value;
            val = obj.getVariable(['m9' num2str(layer) '.' optionValue]);
            BG = obj.convertBG(queriedForm, layer, val);
        end
        
        function BG = convertBG(obj, queriedForm, layer, val)
            l = layer;
            if strcmp(queriedForm,'density')
                % density here is defined as the projected density.
                % this unit is locs/um^2
                optionValue = obj.advanceSetting.(['m9' num2str(l) '_background']).value;
                if strcmp(optionValue,'weight')
                    weight = val;
                    numOfBGLocs = weight*obj.numOfLocsPerLayer(l);
                    BG = numOfBGLocs/(pi*(0.5*obj.roiSize/1000)^2);
                else
                    BG = val;
                end
            elseif strcmp(queriedForm,'weight')
                optionValue = obj.advanceSetting.(['m9' num2str(l) '_background']).value;
                if strcmp(optionValue,'weight')
                    BG = val;
                else
                    density = val;
                    numOfBGLocs = density*(pi*(0.5*obj.roiSize/1000)^2);
                    BG = numOfBGLocs./obj.numOfLocsPerLayer(l);
                end
            else
                error(strjoin(["'" queriedForm "' is not a valid option. Use either 'density', 'weight'"], ''))
            end
        end
        
%         function form = getBGForm(obj, modelnumber)
%             % This function check the form or the background for a layer
%             %
%             % Usage:
%             %   obj.getBGForm(modelnumber)
%             %
%             % Args:
%             %   modelnumber: a integer indicating the layer. E.g., m9X for
%             %   the Xth layer.
%             %
%             % Returns:
%             %   form: a char of either 'density' or 'weight'.
%             
%             obj.getAllParId
%         end
        
        function getLocsInfo(obj)
            %% Get locs counts per layer
            locs = obj.locs;
            allLocsLayers = unique(locs.layer);
            obj.numOfLocsPerLayer = histcounts(locs.layer, 1:max(allLocsLayers)+1);
            obj.representiveLocprec = grpstats(locs.locprecnm, locs.layer, 'median');
            numOfUsedLocs = sum(obj.numOfLocsPerLayer(obj.allModelLayer));
            obj.weightLayer = zeros(size(obj.numOfLocsPerLayer));
            if obj.useCompensation == true
                numOfLocs_maxLayer = max(obj.numOfLocsPerLayer(obj.allModelLayer));
                obj.compensationFactor = numOfLocs_maxLayer./obj.numOfLocsPerLayer;
                obj.weightLayer(obj.allModelLayer) = 1/length(obj.allModelLayer);
            else
                obj.compensationFactor(obj.allModelLayer) = 1;
                obj.weightLayer(obj.allModelLayer) = obj.numOfLocsPerLayer(obj.allModelLayer)/numOfUsedLocs;
            end
        end
        
        function lParSelector(obj, parameterType, parameterForm)
            % Defines the form of parameters to be used
            % It reacts when a certain form of lPar is selected.
            parameterType = strsplit(parameterType, '_');
            targetModel = parameterType{1};
            targetModel = str2num(targetModel(2:end));
            parameterType = parameterType{2};
            switch parameterType
                case 'orientation'
                    % TODO: different rotation form
                case 'background'
                    switch parameterForm
                        case 'density'
                            lPar2rm = ['m' num2str(targetModel) '.offset.weight'];
                        case 'weight'
                            lPar2rm = ['m' num2str(targetModel) '.offset.density'];
                    end
                    if ~isempty(obj.getVariable(lPar2rm))
                        lParArg = defaultOffset(obj, parameterForm);
                        nPar = length(lParArg.name);
                        lParArg.model = repelem(targetModel,nPar);
                        lParArg.type = repelem({'offset'},nPar);
                        lParArg.label = repelem({''},nPar);
                        lParArg.value = obj.getBG(parameterForm,targetModel-90);
                        obj.rmPar(lPar2rm)
                        obj.addPar(lParArg)
                    else
                        warning(['No parameter with the ID: ' lPar2rm '.'])
                        return
                    end
            end
        end
        
        function initLParSelector(obj, modelnumber, isBackground)
            % Initiates the form of parameters to be used
            if isBackground
            	default = defaultLParSelection('background');
                
            else
                default = defaultLParSelection('orientation');
            end
            fn = fieldnames(default);
            for k = 1:length(fn)
                default.(fn{k}).name = ['m' num2str(modelnumber) '_' default.(fn{k}).name];
                default = renameStructField(default,fn{k},['m' num2str(modelnumber) '_' fn{k}]);
                % Check the form of the background: density or weight?
%                 obj.getBGForm(modelnumber);
            end
            obj.addAdvanceSettings(default);
        end

        %% convert related
        function rmConvertRules(obj)
            obj.converterRules.target = {};
            obj.converterRules.rule = {};
        end
        
        function rmOneConvertRule(obj, targets)
            indRm = ismember(obj.converterRules.target,targets);
            obj.converterRules.target(indRm) = [];
            obj.converterRules.rule(indRm) = [];
        end
        
        %% advanced settings
        function initAdvanceSetting(obj)
            default = defaultAdvanceSettings();
            obj.addAdvanceSettings(default);
        end
        
        function updateAdvancedSetting(obj)
            default = defaultAdvanceSettings();
            fn = fieldnames(default);
            for k = 1:length(fn)
                value = obj.getAdvanceSetting(fn{k});
                if isempty(value)
                    % Only add the setting when it doesn't exist
                    obj.addAdvanceSettings(struct(fn{k},default.(fn{k})));
                end
            end
        end
        
        function addAdvanceSettings(obj, settings)
            fn = fieldnames(settings);
            for k = 1:length(fn)
                oneField = settings.(fn{k});
                obj.advanceSetting.(fn{k}) = oneField;
            end
        end
        
        function setAdvanceSetting(obj, par, value)
            obj.advanceSetting.(par).value = value;
        end
        
        function out = getAdvanceSetting(obj, settingName,field)
            try
                if isfield(obj.advanceSetting, settingName)
                switch nargin 
                    case 2
                        out = obj.advanceSetting.(settingName).value;
                    case 3
                        out = obj.advanceSetting.(settingName).(field);
                    otherwise
                        warning("The property 'advanceSetting' might be wrong or missing.")
                        return
                end
                else
                    out = [];
                end
            catch
                warning("The property 'advanceSetting' might be wrong or missing.")
            end
        end
        %% Site registration
        function locs = locsRegister(obj, locs, lParsVal, modelID, varargin)
            if ~isempty(lParsVal)
                locs = obj.locsHandler(locs, lParsVal, modelID, varargin{:});
            end
            % For alignment, get the post-fit transformation
            if ~isempty(obj.converterRules)
                ind = find(startsWith(obj.converterRules.target, 'post_'));
                for k = 1:length(ind)
                    target = obj.converterRules.target{ind(k)};
                    val = obj.getVariable(target);

                    %% Apply the transformation.
                    % Act differenctly depends on the target.
                    allPar = {'x','y','z','scale','zrot'};

                    indPostFit = find(strcmp(strcat('post_', allPar), target));

                    switch indPostFit
                        case 1
                            locs.xnm = locs.xnm+val;
                        case 2
                            locs.ynm = locs.ynm+val;
                        case 3
                            locs.znm = locs.znm+val;
                        case 4
                            locs.xnm = locs.xnm.*val;
                            locs.ynm = locs.ynm.*val;
                            if obj.dataDim == 3
                                locs.znm = locs.znm.*val;
                            end
                        case 5
                            val = deg2rad(val);
                            [locs.xnm,locs.ynm] = rotcoord2(locs.xnm, locs.ynm, val);
                    end
                end

                allPar = {'x','y','z','scale'};
            end
        end
        %% For simlulation
        function assignParsVal(obj)
            % In the given range, assign values to the parameters randomly.
            %
            % Usage:
            %   obj.assignParsVal
            %
            % Args:
            %   obj: an LocMoFit object.
            %
            % Returns:
            %   Parameter values saved in obj.allParsArg.
            lSampling = ~obj.allParsArg.fix;
            rangeSampling = obj.allParsArg.ub(lSampling)-obj.allParsArg.lb(lSampling);
            randVal = rand([sum(lSampling) 1]);
            valueSampling = randVal.*rangeSampling+obj.allParsArg.lb(lSampling);
            obj.allParsArg.value(lSampling) = valueSampling;
        end
        
        function p = getSimIntensity(obj, label)
            obj.assignParsVal;
            parFix = obj.allParsArg.fix;
            obj.allParsArg.fix = true(size(parFix));
            p = obj.intensityCal([], label)';
            obj.allParsArg.fix = parFix;
        end
        
        function modCoord = getSimRef(obj, varargin)
            % Get simulation reference based on the allArgVal.
            %
            % Usage:
            %   modCoord = obj.getSimRef(Name-value)
            %
            % Args:
            %   obj: an LocMoFit object.
            %   Name-value pairs:
            %       * 'finalROISize': the final ROI size for fitting. the 
            %       ROI size for simulations is usually larger than this
            %       value to make sure the background fills everywhere.
            %   
            % Returns:
            %   modCoord: reference coordinates.
            %
            p = inputParser;
            p.addParameter('finalROISize',obj.roiSize)
            p.parse(varargin{:})
            
            finalROISize = p.Results.finalROISize;
                        
            % In the given range, assign values to the parameters randomly.
            obj.assignParsVal;
            
            % Get parameters from convert
            obj.convertNow;
            
            % Get the labels of the model.
            for k = 1:obj.numOfModel
                sigma{k} = obj.model{k}.sigma;
                sigmaSet{k} = obj.model{k}.sigmaSet;
                obj.model{k}.sigmaSet = [];
                obj.model{k}.sigma = 5;
            end
            
            %% Get the origainl model
            modCoord = obj.getLayerPoint(0.75);
            
            for k = 1:obj.numOfModel
                obj.model{k}.sigma = sigma{k};
                obj.model{k}.sigmaSet = sigmaSet{k};
            end
            
            modCoord = transModPoint(obj, modCoord);
            
            for l = 1:length(obj.allModelLayer)                
                if any(modCoord{l}.n>1)
                    % locs.n scales to the number of labels at
                    % corresponding positions
                    repeatLabel = round(modCoord{l}.n./min(modCoord{l}.n));
                    labelID = 1:length(modCoord{l}.n);
                    allLabels = repelem(labelID, repeatLabel);
                    numOfLabels = length(allLabels);
                    numOfMol = obj.getVariable(['m' num2str(90+obj.allModelLayer(l)) '.numOfMol']);
                    if ~isempty(numOfMol)&&numOfMol>0
%                         numOfLabels = numOfMol;
                        indKept = randi(numOfLabels, numOfMol,1);
                        label_kept = allLabels(indKept);
                        modCoord{l}.x = modCoord{l}.x(label_kept);
                        modCoord{l}.y = modCoord{l}.y(label_kept);
                        modCoord{l}.n = modCoord{l}.n(label_kept);
                        if obj.dimension == 3
                            modCoord{l}.z = modCoord{l}.z(label_kept);
                        end
                        numOfLabels = numOfMol;
                    end
                else
                    numOfLabels = length(modCoord{l}.x);
                    numOfMol = obj.getVariable(['m' num2str(90+obj.allModelLayer(l)) '.numOfMol']);
                    if ~isempty(numOfMol)&&numOfMol>0
                        numOfLabels = numOfMol;
                        indKept = randi(length(modCoord{l}.x), numOfMol,1);
                        modCoord{l}.x = modCoord{l}.x(indKept);
                        modCoord{l}.y = modCoord{l}.y(indKept);
                        modCoord{l}.n = modCoord{l}.n(indKept);
                        if obj.dimension == 3
                            modCoord{l}.z = modCoord{l}.z(indKept);
                        end
                    end
                end
                
                
                %% Assign background for each layer
                % 1. First check the number of labels in the model.
                % 2. And then based on the weight of offset (background),
                % calculate the probability of getting a label in the
                % background.
                % 3. Randomly assign the x y z coordinates to the offset labels
                
                
                
                % show the diff between exbect and experimental
                %             round((numOfLabels*pExp_Offset)/(1-pExp_Offset))-numOfLabels_offset
                
                if 0
                    % circular ROI
                    % 210102: not used anymore
                    
                    pExp_Offset = obj.getVariable(['par.m9' num2str(l) '.offset.weight']);
                    numAll = round(numOfLabels*(pExp_Offset)/(1-pExp_Offset)+numOfLabels);
                    pObs_Offset = rand([numAll 1]);
                    lKept = pExp_Offset>pObs_Offset;
                    numOfLabels_offset = sum(lKept);
                    
                    u = rand([numOfLabels_offset 1])+rand([numOfLabels_offset 1]);
                    u(u>1) = 2-u(u>1);
                    r = u;
                    theta = rand([numOfLabels_offset 1]) * 2 * pi;
                    
                    x_offset = r .* cos(theta);
                    y_offset = r .* sin(theta);
                    
                    roiSize = obj.roiSize;
                    
                    modCoord{l}.x = [modCoord{l}.x; x_offset.*roiSize/2];
                    modCoord{l}.y = [modCoord{l}.y; y_offset.*roiSize/2];
                    n = ones(size(x));
                    modCoord{l}.n = [modCoord{l}.n; n];
                else
                    % cilinderical ROI in the cubic
                    % scale up to ROI size of the simulation
                    
                    numOfLabels;
                    
                    offset = obj.exportOffset();
                    if isfield(offset{90+l},'density')
                        BGDensity = offset{90+l}.density;
                        BGCount = BGDensity*(pi*(0.5*finalROISize/1000)^2);
                        offset{90+l}.weight = BGCount/(BGCount+numOfLabels);
                    end
                    pExp_Offset = offset{90+l}.weight;
                    
                    numAll = numOfLabels*(pExp_Offset)/(1-pExp_Offset)+numOfLabels;
                    
                    roiSize = obj.roiSize;
                    [x,y,z] = getBGLocs(numAll, roiSize, finalROISize, pExp_Offset, obj.dataDim);
                    n = ones(size(x));
                    modCoord{l}.x = [modCoord{l}.x; x];
                    modCoord{l}.y = [modCoord{l}.y; y];
                    modCoord{l}.n = [modCoord{l}.n; n];
                end
                if ~isempty(z)
                    % for 3D model
                    modCoord{l}.z = [modCoord{l}.z; z];
                end
            end
        end
        function distribution = getLLExpDist(obj, n)
            % Get the labels of the model.
            if ~strcmp(obj.model{1}.modelType,'image')
                for k = 1:obj.numOfModel
                    sigma{k} = obj.model{k}.sigma;
                    sigmaSet{k} = obj.model{k}.sigmaSet;
                    obj.model{k}.sigmaSet = [];
                    obj.model{k}.sigma = 5;
                end
                
                %% Get the origainl model
                modCoord = obj.getLayerPoint(0.75);
                modCoord = transModPoint(obj, modCoord);
                
                for k = 1:obj.numOfModel
                    obj.model{k}.sigma = sigma{k};
                    obj.model{k}.sigmaSet = sigmaSet{k};
                end
                
                for l = 1:obj.numOfModel
                    nLocs = sum(obj.locs.layer==l);
                    bgWeight = obj.exportPars(91,'offset');
                    var = obj.getVariable('m1.lPar.variation');
                    nModCoord = length(modCoord{l}.x);
                    for nn = n:-1:1
                        propNonBG = rand([nLocs 1]);
                        lNonBG = propNonBG>bgWeight.weight;
                        nLoc_model = sum(lNonBG);
                        locprecnm = obj.locs.locprecnm(lNonBG);
                        locprecznm = obj.locs.locprecznm(lNonBG);
                        locprecnm_ = sqrt(locprecnm.^2+var^2);
                        locprecznm_ = sqrt(locprecznm.^2+var^2);
                        
                        idxModCoord = randi(nModCoord, [nLoc_model 1]);
                        rsModCoord.xnm = modCoord{l}.x(idxModCoord)+randn([nLoc_model 1]).*locprecnm_;
                        rsModCoord.ynm = modCoord{l}.y(idxModCoord)+randn([nLoc_model 1]).*locprecnm_;
                        rsModCoord.znm = modCoord{l}.z(idxModCoord)+randn([nLoc_model 1]).*locprecznm_;
%                         rsModCoord.xnm = modCoord{l}.x(idxModCoord);
%                         rsModCoord.ynm = modCoord{l}.y(idxModCoord);
%                         rsModCoord.znm = modCoord{l}.z(idxModCoord);
                        rsModCoord.locprecnm = locprecnm;
                        rsModCoord.locprecznm = locprecznm;
                        
                        if sum(~lNonBG)>0
                            locprecnm_BG = obj.locs.locprecnm(~lNonBG);
                            locprecznm_BG = obj.locs.locprecznm(~lNonBG);
                            
                            roiSize = obj.roiSize;
                            finalROISize = obj.roiSize;
                            [x,y,z] = getBGLocs(length(locprecnm_BG), roiSize, finalROISize, 1, obj.dataDim);
                            rsModCoord.xnm = [rsModCoord.xnm; x];
                            rsModCoord.ynm = [rsModCoord.ynm; y];
                            
                            if ~isempty(z)
                                rsModCoord.znm = [rsModCoord.znm; z];
                            end
                            rsModCoord.locprecnm = [rsModCoord.locprecnm; locprecnm_BG];
                            rsModCoord.locprecznm = [rsModCoord.locprecznm; locprecznm_BG];
                        end
                        
                        rsModCoord.layer = ones([nLocs 1]);
                        distribution(nn) = obj.loglikelihoodFun([],[],rsModCoord)/nLocs;
                    end
                    
                end
            else
                distribution = [];
            end
        end
        function setModelSetting(obj, modelInd, option, value)
            obj.model{modelInd}.(option) = value;
        end
        function updateFitInfo(obj)
            [obj.fitInfo.LLfit, obj.fitInfo.weightModel] = intensityCal(obj, fitPars, locs);
        end
        %% Temporary info.
        function setTemp(obj, name, val)
            obj.temp.(name) = val;
        end
        function val = getTemp(obj, name)
            if isfield(obj.temp, name)
                val = obj.temp.(name);
            else
                val = [];
            end
        end
        function rmTemp(obj, name)
            obj.temp = rmfield(obj.temp, name);
        end
        
        %% helper function
        function pars = fitTest(obj, parID)
            oRefPoint_spacing = obj.refPoint_spacing;
            oAllParsArg = obj.allParsArg;
            samplingSpacing = [1.5 1.25 1 0.75 0.5 0.25 0.2];
            for l = length(samplingSpacing):-1:1
                obj.refPoint_spacing = samplingSpacing(l);
                obj.fit(obj.locs);
                for k = length(parID):-1:1 
                    pars{k}(l) = obj.getVariable(parID{k});
                end
            end
            obj.refPoint_spacing = oRefPoint_spacing;
            obj.allParsArg = oAllParsArg;
        end
        
        function saveHandles(obj, pard, tag)
%             if isempty(varargin) || rem(length(varargin),2)>0
%                 error('Wrong pair(s) of field names and handles.')
%             end
%             fn = varargin(1:2:end);
%             h = varargin(2:2:end);
%             for k = 1:length(fn)
%                 obj.handles.(fn{k}) = h{k};
%             end
            fn = fieldnames(pard);
            for k = 1:length(fn)
                if exist('tag', 'var')
                    obj.handles.([fn{k} '_' tag]) = pard.(fn{k});
                else
                    obj.handles.(fn{k}) = pard.(fn{k});
                end
            end
        end
        
        function rmHandles(obj, fn)
            obj.handles = rmfield(obj.handles,fn);
        end
        
        function LLfit = loglikelihoodFun(obj, fitPars, compensationFactor, locs, varargin)
            % loglikelihoodFun computes the log-likelihood value of the
            % fit.
            obj.currentFitPars = fitPars;
            p = inputParser;
            p.addParameter('expected',false);
            p.addParameter('dx',1);
            p.parse(varargin{:});
            results = p.Results;
            
            if ~exist('fitPars', 'var')||isempty(fitPars)
                fitPars = obj.allParsArg.value(~obj.allParsArg.fix)';
            end
            
            if ~exist('locs', 'var')||isempty(locs)
                locs = obj.locs;
            end
            
            if ~exist('compensationFactor', 'var')||isempty(compensationFactor)
                % Get compensationFactor
                % This factor compensate the number of locs between different channels
                compensationFactor = zeros(size(locs.layer))';
                for k = 1:obj.numOfLayer
                    compensationFactor(locs.layer==k) = obj.compensationFactor(k);
                end
                compensationFactor(~ismember(locs.layer, obj.allModelLayer)) = [];
            end
            
            if results.expected    
                % expected LL
                gridPos = locs;
                sizeFields = size(gridPos.xnm);
                for k = length(obj.locs.xnm):-1:1
                    gridPos.locprecnm = repmat(obj.locs.locprecnm(k),sizeFields);
                    gridPos.locprecznm = repmat(obj.locs.locprecznm(k),sizeFields);
                    prob = obj.intensityCal(fitPars,gridPos);
                    % go through all the locs
                    LLfit(k) = sum(prob.*log(prob),2)/sum(prob);
                end
                LLfit = mean(LLfit);
                if 0
                    prob = obj.intensityCal(fitPars,locs);
                    LLfit = sum(prob.*log(prob),2)/sum(prob);
                end
            else
                prob = obj.intensityCal(fitPars,locs);
                if isempty(prob)
                    'hey'
                end
                LLfit = sum(compensationFactor.*log(prob),2);
            end
            if strcmp(obj.status, 'initial')
                obj.status = 'iterative';
            end
        end
        
        function ELL = getELL(obj, fitPars, compensationFactor, dx, varargin)
            p = inputParser;
            p.addParameter('representiveLocprec',obj.representiveLocprec);
            p.addParameter('representiveLocprecz',obj.representiveLocprec*3);
            p.parse(varargin{:});
            results = p.Results;
            
            gridPos1D = -obj.roiSize/2:dx:obj.roiSize/2;
            
            switch obj.dataDim
                case 2
                    [X, Y] = meshgrid(gridPos1D,gridPos1D);
                    gridPos.xnm = X(:);
                    gridPos.ynm = Y(:);
                case 3
                    [X, Y, Z] = meshgrid(gridPos1D,gridPos1D,gridPos1D);
                    gridPos.xnm = X(:);
                    gridPos.ynm = Y(:);
                    gridPos.znm = Z(:);
            end
            
            lOutOfRange = sqrt(gridPos.xnm.^2+gridPos.ynm.^2)>obj.roiSize/2;
            gridPos.xnm(lOutOfRange)=[];
            gridPos.ynm(lOutOfRange)=[];
            gridPos.znm(lOutOfRange)=[];
            
            originalSize = size(gridPos.xnm);
            l = obj.numOfLayer;
            if l > 1
                gridPos.xnm = repmat(gridPos.xnm, l, 1);
                gridPos.ynm = repmat(gridPos.ynm, l, 1);
                if obj.dataDim == 3
                    gridPos.znm = repmat(gridPos.znm, l, 1);
                end
            end
            gridPos.layer = zeros(originalSize);
            gridPos.layer = repmat(gridPos.layer, 1, l)+(1:l);
            gridPos.layer = gridPos.layer(:);
            
            gridPos.locprecnm = repmat(results.representiveLocprec,originalSize);
            gridPos.locprecznm = repmat(results.representiveLocprecz,originalSize);
            ELL = obj.loglikelihoodFun(fitPars, compensationFactor, gridPos, 'expected', true, 'dx', dx);
        end
        
        function OFLL = getOFLL(obj,compensationFactor)
            % get overfitted log-likelihood
            locs = obj.locs;
            model = locsModel(locs);
            model.dimension = 3;
            ctrlFitter = LocMoFit;
            ctrlFitter.dataDim = 3;
            ctrlFitter.addModel(model);
            ctrlFitter.setParArg('m1.lPar.x','fix',true);
            ctrlFitter.setParArg('m1.lPar.y','fix',true);
            ctrlFitter.setParArg('m1.lPar.z','fix',true);
            ctrlFitter.setParArg('m1.lPar.zrot','fix',true);
            ctrlFitter.setParArg('m1.lPar.xrot','fix',true);
            ctrlFitter.setParArg('m1.lPar.yrot','fix',true);
            ctrlFitter.setParArg('m1.lPar.weight','fix',true);
            ctrlFitter.setParArg('m91.offset.weight','fix',true);
            
            OFLL = ctrlFitter.loglikelihoodFun([],compensationFactor,locs)/size(locs.xnm,1);
        end
        
        function out = rel(obj,value,posIn,posOut)
            xrot = obj.getVariable('pars.m1.lPar.xrot');
            yrot = obj.getVariable('pars.m1.lPar.yrot');
            zrot = obj.getVariable('pars.m1.lPar.zrot');
            if obj.dataDim == 3
                pos = zeros([1 3]);
                pos(posIn) = value;
                [out(1),out(2),out(3)] = rotcoord3(pos(1),pos(2),pos(3),deg2rad(xrot),deg2rad(yrot),deg2rad(zrot),'XYZ');
            else
                pos = zeros([1 2]);
                disp('The function rel for 2D needs to be implemented.')%!!!!
            end
            out = out(posOut);
        end
        
        %% Geometric model related.
        function derivedPars = getDerivedPars(obj, varargin)
            % Get derived parameters of all (default) or a specific model.
            % Args:
            % 	obj: an :class:`LocMoFit` object.
            % Returns:
            %   
            %
            if isempty(varargin)
                for m = obj.numOfModel:-1:1
                    allPars = obj.exportPars(m, 'allPar');
                    if isempty(obj.model{m}.modelObj)
                        derivedPars{m} = [];
                    else
                        derivedPars{m} = obj.model{m}.modelObj.getDerivedPars(allPars);
                    end
                end
            end
            if isempty(derivedPars)
                obj.fitInfo.derivedPars = [];
            else
                obj.fitInfo.derivedPars = derivedPars;
            end
            
        end
        
        function settings = getModelInternalSettingList(obj, modelID)
            % Get the list of the internal settings' name.
            %
            % --- Syntax ---
            % settings = getModelInternalSettingList(obj, modelInd).
            % 
            % --- Description---
            % settings: a list (string array) of the internal settings' name.
            % obj: an LocMoFit object.
            % modelID: the ID of the model.
            if ~isempty(obj.model{modelID}.modelObj)&&~isempty(obj.model{modelID}.modelObj.internalSettings)
                settings = fieldnames(obj.model{modelID}.modelObj.internalSettings);
            else
                settings = [];
            end
        end
        function setModelInternalSetting(obj, modelInd, setting, value)
            % Set the value of a model internal settings.
            if ~isempty(obj.model{modelInd}.modelObj)
                obj.model{modelInd}.modelObj.setInternalSettings(setting,value);
            end
        end
        function value = getModelInternalSetting(obj, modelInd, setting)
            % Get the value of a model internal settings.
            if ~isempty(obj.model{modelInd}.modelObj)
                value = obj.model{modelInd}.modelObj.getInternalSettings(setting);
            else
                value = [];
            end
        end
        function mParsArgModified_callback(obj, b)
             obj.updateModel(b.modelID)
        end
        
        %% get/set methods for parameters
        function val = get.numOfModel(obj)
            val = length(obj.model);
        end
        
        function val = get.numOfLayer(obj)
            val = length(obj.allModelLayer);
        end
        
        function val = get.roiAreaVol(obj)
            switch obj.dataDim
                case 2
                    val = pi*obj.roiSize^2;
                case 3
                    % cylindrical
                    % TODO: better way for estimating the height
                    val = pi*obj.roiSize^2*obj.roiSize;
            end
        end
        
        function set.advanceSetting(obj, val)
            oldVal = obj.advanceSetting;
            fn = fieldnames(val);
            % detect model specific settings
            indModelSpecific = find(~cellfun(@isempty, (regexp(fn,'^(m\d+\_)','start'))));
            if ~isempty(indModelSpecific)
                for k = 1:length(indModelSpecific)
                    oneFn = fn{indModelSpecific(k)};
                    if isfield(oldVal, oneFn)&&~isequal(oldVal.(oneFn).value,val.(oneFn).value)
                    	obj.lParSelector(oneFn, val.(oneFn).value)
                    end
                end
            end
            obj.advanceSetting = val;
        end
        %% for compatibilities
        function updateVersion(obj)
            % This function is for necessary updates
            % Structral changes leading to running failular have to be
            % fixed here.
            uModelTypeOptoins(obj)
            uEarlier(obj)
            u210630(obj)
            obj.updateAdvancedSetting
%             u210917(obj)
%             u210919(obj)            
        end
    end
    methods(Access = protected)
      % Override copyElement method:
      function cp = copyElement(obj)
         cp = copyElement@matlab.mixin.Copyable(obj);
         for m = 1:cp.numOfModel
            cp.model{m} = copy(cp.model{m});
            cp.model{m}.ParentObject = cp;
         end
      end
    end
    methods (Static)
        function xPars = pars2struct(fn, xPars, sourcePars, k, varargin)
            p = inputParser;
            p.addParameter('AddUp',1);
            parse(p, varargin{:})
            results = p.Results;
            addUp = results.AddUp;
            %%
            % xPars = pars2struct(fn, xPars, sourcePars, k)
            % xPars and sourcePars are numeric array.
            % For kth model, sum up the value of xPars and sourcePars, and save
            % them to the fn fields of the struct xPars
            for l = 1:length(fn)
                if k>1 & addUp == 1
                    xPars{k}.(fn{l}) = xPars{1}.(fn{l}) + sourcePars(:,l)';     % lPars are all relative to the model one
                else
                    xPars{k}.(fn{l}) = sourcePars(:,l)';
                end
            end
        end
    end
    events
        mParsArgModified
    end
end

function lPars = defaultLpars(obj, dim)
    % get default lPars
    % define the default lPars here

    lPars.name = {'x' 'y' 'zrot' 'variation' 'xscale' 'yscale' 'weight'}';
    lPars.fix = logical([0 0 0 1 1 1 0]');
    lPars.lb = [-obj.roiSize/2 -obj.roiSize/2 -inf 0 1 1 1e-5]';
    lPars.ub = [obj.roiSize/2 obj.roiSize/2 inf 10 1 1 1]';
    lPars.min = [-obj.roiSize/2 -obj.roiSize/2 -inf 0 1 1 1e-5]';
    lPars.max = [obj.roiSize/2 obj.roiSize/2 inf 20 1 1 1]';
    lPars.value = [0 0 0 0 1 1 1e-5]';
    if dim == 3
        lPars.name = [lPars.name; {'z' 'xrot' 'yrot' 'zscale'}'];
        lPars.fix = [lPars.fix; logical([0 0 0 1])'];
        lPars.lb = [lPars.lb; [-obj.roiSize -inf -inf 1]'];
        lPars.ub = [lPars.ub; [obj.roiSize inf inf 1]'];
        lPars.min = [lPars.min; [-obj.roiSize -inf -inf 1]'];
        lPars.max = [lPars.max; [obj.roiSize inf inf 1]'];
        lPars.value = [lPars.value; [0 0 0 1]'];
    end
end

function offset = defaultOffset(obj, name)% get default lPars
    % define the default lPars here

    offset.name = {'weight', 'density'}';
    offset.fix = [false false];
    offset.lb = [0.001 1];
    offset.ub = [0.999 1000];
    offset.min = [0 0];
    offset.max = [1 inf];
    offset.value = [0 0];
   
    lExoprt = ismember(offset.name, name);
    fn = fieldnames(offset);
    for k = 1:length(fn)
        offset.(fn{k}) = offset.(fn{k})(lExoprt);
    end
end

function out = defaultAdvanceSettings
    % This function defines the default values of advance settings.
    % The function initAdvanceSettings(obj) calls this function.
    out.controlLogLikelihood.option = {'none','expected','overfitted','both'};
    out.controlLogLikelihood.value = 'none';
    out.controlLogLikelihood.name = 'Control log-likelihood';

    out.gaussDistMode.option = {'ordinary', 'fast'};
    out.gaussDistMode.value = 'ordinary';
    out.gaussDistMode.name = 'Gauss distance';
    
    out.gaussDistCutoff.option = {};
    out.gaussDistCutoff.value = 3.5;
    out.gaussDistCutoff.name = 'Gauss distance offset';
    
    out.confidenceInterval.option = {'on','off'};
    out.confidenceInterval.value = 'off';
    out.confidenceInterval.name = 'Confidence interval';
end

function [out,allOptions] = defaultLParSelection(parameterType)
    % This function defines the default values of the lPar selections.
    % The function initLParSelector(obj) calls this function.
    if exist('parameterType','var')
        switch parameterType
            case 'orientation'
                out.orientation.option = {'ZYX','zxz'};
                out.orientation.value = 'ZYX';
                out.orientation.name = 'Orientation';
            case 'background'
                out.background.option = {'weight', 'density'};
                out.background.value = 'weight';
                out.background.name = 'Background';
            otherwise
                warning(['defaultLParSelection: The parameter type "' parameterType '" is not supported.'])
                return
        end
    else
        out = [];
    end
    if nargout>1
        allOptions = {'orientation','background'};
    end
end

%% For version check and update
function u210917(obj)
    fn = fieldnames(obj.advanceSetting);
    numOfModel = obj.numOfModel;
    allOptions = defaultXXX(obj, {'weight','density'}); % yet to be implemented
    allOptionsName = allOptions.name;
    if 0
        for m = 1:numOfModel
            allOptionsName_m = strcat('m', num2str(m), '_', allOptionsName);
            if any(ismember(allOptionsName_m,fn))
                obj.initLParSelector(m, false)
            end
        end
    end
end

function u210919(obj)
    if isempty(obj.modelLayer)
        obj.modelLayer = getFieldAsVector(obj.model,'layer');
    end
end

function u210630(obj)
    fn = fieldnames(obj.advanceSetting);
    layers = unique(obj.allModelLayer);
    [~,allOptions] = defaultLParSelection;
    
    for l = layers
        allOptions_l = strcat('m', num2str(90+l), '_', allOptions);
        if ~any(ismember(allOptions_l,fn))
            obj.initLParSelector(90+l, true)
        end
    end
end

function uEarlier(obj)
    if isempty(obj.advanceSetting)
        obj.initAdvanceSetting;
    end
end

function uModelTypeOptoins(obj)
    for m = 1:obj.numOfModel
        oneModelObj = obj.model{m}.modelObj;
        if ~isempty(oneModelObj)
            modelTypeOption = oneModelObj.modelTypeOption;
            if isempty(modelTypeOption)
                modClass = class(oneModelObj);
                modObjContainer = eval(modClass);
                defaultModelTypeOption = modObjContainer.modelTypeOption;
                oneModelObj.modelTypeOption = defaultModelTypeOption;
                currentModType = oneModelObj.modelType;
                if ~strcmp(currentModType, defaultModelTypeOption)
                    if strcmp(currentModType,'discrete')&&any(strcmp('discretized', defaultModelTypeOption))
                        oneModelObj.modelType = 'discretized';
                        oldModType = obj.model{m}.modelType;
                        obj.model{m}.modelType = 'discretized';
                        warning(['Model ' num2str(m) ': its model type was changed from [' oldModType '] to [discretized] due to an update. If this is not expected, check the model type carefully.'])
                    else
                        warning(['Model ' num2str(m) ': ambiguous model type. Check whether the model type is outdated or not and consider to update it.'])
                    end
                end
                oneModelObj.modelTypeOption = modObjContainer.modelTypeOption;
            end
        end
    end
end
% 
% 
% hold(subax1, 'on')
% axes(subax1);
% plot(subax1, (newlocs.xnmrot+150),(newlocs.ynmrot+150), ' xk', 'MarkerSize', 10, 'LineWidth', 2)
% hold(subax2, 'on')
% plot(subax2, (newlocs.znm+150),(newlocs.xnmrot+150), ' xk', 'MarkerSize', 10, 'LineWidth', 2)