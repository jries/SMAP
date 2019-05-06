classdef SMLMModelFitter<interfaces.SEEvaluationProcessor
    % This is a plugin in development. Public has no access to the 
    % run_ functions called in this plugin. For internal user, you need
    % NPC3D to run this plugin.
    methods
        function obj=SMLMModelFitter(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        
        function out=run(obj,p)
            out = runSMLMModelFitter(obj,p);
        end
        function pard=guidef(obj)
            %init
            obj.setPar('userDefinedPars', userDefinedParsTemp);
            addpath('..\smlm_geomodel_fitter', '..\smlm_geomodel_fitter\external')
            pard=guidef(obj);
        end
    end
end

function pard = guidef(obj)
    pard.t1.object = struct('Style','text','String', 'Number of Structures:');
    pard.t1.position = [1 1];
    pard.t1.Width = 2;
    
    pard.numOfStructures.object = struct('Style','edit','String','1');
    pard.numOfStructures.position = [1 3];
    pard.numOfStructures.Width = 1.5;
    
    pard.setNumOfStruct.object = struct('Style','pushbutton','String','Set', 'Callback', {{@setNumOfStruct_callback,obj}});
    pard.setNumOfStruct.position = [1 4.5];
    pard.setNumOfStruct.Width = 0.5;
    
    pard.t2.object = struct('Style','text','String', 'Current step');
    pard.t2.position = [2 1];
    pard.t2.Width = 1;
    
    pard.currentStep.object = struct('Style','popupmenu','String', {'1'}, 'Value', 1);
    pard.currentStep.position = [2 2];
    pard.currentStep.Width = 0.5;
    
    pard.setModel.object = struct('Style','pushbutton','String', 'Select models', 'Callback', {{@selectModels_callback,obj}});
    pard.setModel.position = [3 1];
    pard.setModel.Width = 1;
    
    pard.setIPar.object = struct('Style','pushbutton','String', 'Intrinsic Par', 'Callback', {{@setIPar_callback,obj}});
    pard.setIPar.position = [3 2];
    pard.setIPar.Width = 1;
    
    pard.setEPar.object = struct('Style','pushbutton','String', 'Extrinsic Par', 'Callback', {{@setEPar_callback,obj}});
    pard.setEPar.position = [3 3];
    pard.setEPar.Width = 1;
    
    pard.reset.object = struct('Style','pushbutton','String', 'Reset', 'Callback', {{@reset_callback,obj}});
    pard.reset.position = [3 4];
    pard.reset.Width = 1;
    
    pard.displayFolder.object = struct('Style','edit','String', '---', 'Enable','off');
    pard.displayFolder.position = [4 1];
    pard.displayFolder.Width = 3;
    
    pard.modelFolder.object = struct('Style','pushbutton','String', 'Model', 'Callback', {{@load_callback,obj}});
    pard.modelFolder.position = [4 4];
    pard.modelFolder.Width = 1;
    
    pard.symbol.object = struct('Style','edit','String', 'a');
    pard.symbol.position = [5 1];
    pard.symbol.Width = 0.5;
    
    pard.linkPar.object = struct('Style','popupmenu','String', {'none'}, 'Value', 1, 'Enable', 'off');
    pard.linkPar.position = [5 1.5];
    pard.linkPar.Width = 2;
    
    pard.assign.object = struct('Style','pushbutton','String', 'Assign', 'Callback', {{@assign_callback,obj}});
    pard.assign.position = [5 3.5];
    pard.assign.Width = 1;
    
    pard.symbol2.object = struct('Style','edit','String', 'e.g., meanAB');
    pard.symbol2.position = [6 1];
    pard.symbol2.Width = 0.5;
    
    pard.expression.object = struct('Style','edit','String', 'e.g., mean([a b])');
    pard.expression.position = [6 1.5];
    pard.expression.Width = 2;
    
    pard.defineNewPar.object = struct('Style','pushbutton','String', 'Define', 'Callback', {{@defineNewPar_callback,obj}});
    pard.defineNewPar.position = [6 3.5];
    pard.defineNewPar.Width = 1;
    
    pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi'};
    pard.plugininfo.type='ROI_Evaluate';
end

function setIPar_callback(a,b,obj)
    % get basis
    [~, currentStep, numOfStep, numOfStructure] = getParInfo(obj);
    
    % GUI part
    % one structure one tab
    [~,~,tab] = createParGui(currentStep,numOfStructure,'Intrinsic parameters');
    
    % data part
    defaultPar(obj, 'lPar', numOfStep, numOfStructure, intrinsicParTemp())
    lPar = obj.getPar('lPar');
	oneIPar = lPar{currentStep};
    
    allUit = {};
    for j = 1:numOfStructure(currentStep)
        uit = uitable(tab{j},'Data',oneIPar{j}, 'Position', [10 10 520 360]);
        uit.ColumnEditable = [false true true true true true true];
        uit.CellEditCallback = {@updatePars,j};
        allUit{j} = uit;
    end
    
    function updatePars(src,event,j)
        lPar{currentStep}{j} = allUit{j}.Data;
        obj.setPar('lPar',lPar)
    end
end

function setEPar_callback(a,b, obj)
    % get basis
    [~, currentStep, numOfStep, numOfStructure] = getParInfo(obj);
    
    % GUI part
    % one structure one tab
    [~,~,tab] = createParGui(currentStep, numOfStructure, 'Extrinsic parameters');
    allModels = obj.getPar('allModels');
    selectedModelIdx = obj.getPar('selectedModel');
    selectedModel = cellfun(@(x)allModels(x,:), selectedModelIdx, 'UniformOutput', false);
%     str2fun(selectedModel)
%     NPCPointModel([],[])
    defaultPar(obj, 'mPar', numOfStep, numOfStructure, selectedModel)
   
    mPar = obj.getPar('mPar');
	oneEPar = mPar{currentStep};

    for j = 1:numOfStructure(currentStep)
        uit = uitable(tab{j},'Data',oneEPar{j}, 'Position',[10 10 520 360]);
        uit.ColumnEditable = [false true true true true true true];
        uit.CellEditCallback = {@updatePars,j};
        allUit{j} = uit;
    end
    
    function updatePars(src,event,j)
        mPar{currentStep}{j} = allUit{j}.Data;
        obj.setPar('mPar',mPar)
    end
end

function [selMol, currentStep, numOfStep, numOfStructure] = getParInfo(obj)
    selMol = obj.getPar('selectedModel');
    currentStep = obj.getSingleGuiParameter('currentStep').Value;
    numOfStep = length(selMol);
    numOfStructure = obj.getSingleGuiParameter('numOfStructures');
end

function [f, tabgp, tab]=createParGui(currentStep, numOfStructure, winTitle)
%     f = figure('Name', winTitle);
    f = uifigure('Name', winTitle);
%     tabgp = uitabgroup(f,'Position',[.01 .01 0.98 0.98]);
    tabgp = uitabgroup(f,'Position',[10 10 540 400]);
    numOfStructureOneStep = numOfStructure(currentStep);
    tab=cell(1:numOfStructureOneStep);
    for k = 1:numOfStructureOneStep
        tab{k} = uitab(tabgp,'Title',num2str(k));
    end
end

function defaultPar(obj, onePar, numOfStep, numOfStructure,parStruct)
    allPars = obj.getPar('allPars');
    if isempty(obj.getPar(onePar))
        data=cell(1,numOfStep);
        for l = 1:numOfStep
            for m = 1:numOfStructure(l)
                % initiate parameter table
                if isstruct(parStruct)
                    pars = fieldnames(parStruct);
                else
                    % for images
                    if contains(parStruct{l}(m,:),'.mat')
                        pars = [];
                    else
                        modelfunc = str2func(parStruct{l}(m,:));
                        pars = modelfunc();
                        pars = fieldnames(pars.par);
                    end
                end
                logicalVal = false(size(pars));
                numericVal = zeros(size(pars));
                
                % mount user defined pars
                userDefinedPars = obj.getPar('userDefinedPars');
                if l>1
                    allPrePars = allPars(allPars.Step<l,:);
                    predefinedOption = [repmat('step', size(allPrePars,1),1) num2str(allPrePars.Step) repmat('.structure', size(allPrePars,1),1) num2str(allPrePars.Structure) repmat('.', size(allPrePars,1),1) char(allPrePars.ParType) repmat('.', size(allPrePars,1),1) char(allPrePars.Parameters)];
                    
                    predefined = categorical(repmat({'none'},size(pars,1), 1),['none'; cellstr(predefinedOption); userDefinedPars.symbol]);% mount predefined and user-defined pars
                else
                    predefined = categorical(repmat({'none'},size(pars,1), 1),['none'; userDefinedPars.symbol]);
                end
                % still need to be checked
                oneData = table(pars, logicalVal, predefined, numericVal, numericVal, numericVal, numericVal);
                oneData.Properties.VariableNames = {'Parameters','Fix','Fix_fromFit','Fix_value','Lower_bound','Upper_bound','Starting'};
                data{l}{m}=oneData;
            end
        end
        obj.setPar(onePar, data);
    end
end

function load_callback(a,b,obj)
    default = obj.getSingleGuiParameter('displayFolder');
    if isequal(default, '---')
        f = './';
    else
        f = default;
    end
    f=uigetdir(f,'Choose the folder of models');
    if ~f
        return
    end
    addpath(f);
    allModels = ls([f '\*.m']);
    allModels = [string(allModels); string(ls([f '\*.mat']))];
    allModels = deblank(allModels);
    allModels = regexprep(allModels,'\.m$','');
    obj.setGuiParameters(struct('displayFolder',f));
    obj.setPar('allModels',allModels);
end

function setNumOfStruct_callback(a,b,obj)
    numOfStructures = obj.getSingleGuiParameter('numOfStructures');
    obj.guihandles.currentStep.String = cellfun(@(x)num2str(x),num2cell(1:length(numOfStructures)), 'UniformOutput', false);
    obj.setPar('selectedModel',[]);
end

function reset_callback(a,b,obj)
    obj.setPar('lPar',[])
    obj.setPar('mPar',[])
    obj.setPar('userDefinedPars',userDefinedParsTemp)
end

function selectModels_callback(a,b,obj)
% the data structure is selMod{l}(m), where l specify the lth step, and the
% m specify the mth structure of the lth step
    f = figure('Name', 'Select Models', 'CloseRequestFcn', {@close_callback, obj});
    selMod = obj.getPar('selectedModel');
    c_step = {};
    c_structure = {};
    if isempty(selMod)
        selMod = {};
    else
    end
    numOfStructure = obj.getSingleGuiParameter('numOfStructures');
    numOfStep = length(numOfStructure);
    for l = 1:numOfStep
        for m = 1:numOfStructure(l)
            if length(selMod)<l
                selMod{l} = [];
            end
            if isempty(selMod)|length(selMod{l})~=numOfStructure(l)
                selMod{l}(m) = 1;
            end
            c_structure{m} = uicontrol(f, 'Style','popupmenu','String',obj.getPar('allModels'),'Value',selMod{l}(m));
            c_structure{m}.Position = [20+(l-1)*160 350-(m-1)*20 150 20];
            c_structure{m}.Callback = {@updatePars,{obj,l,m}};
        end
        c_step{l} = c_structure;
    end
    
    function updatePars(a,b,obj)
        l = obj{2};
        m = obj{3};
        selMod{l}(m) = b.Source.Value;
        obj{1}.setPar('selectedModel',selMod);
    end

    function close_callback(a,b,obj)
        allModels = obj.getPar('allModels');
        numOfStructures = obj.getSingleGuiParameter('numOfStructures');
        numSteps = length(numOfStructures);
        
        % create the par table template
        allPars = table('Size',[0 6],'VariableTypes',{'double','double','string','string','double','string'},'VariableNames',{'Step','Structure','ParType','Parameters','Val','Symbol'});
        locsPar = fieldnames(intrinsicParTemp());
        numLPar= size(locsPar,1);
        mTypes = cell([1 numSteps]);
        
        % pre-locate the cell for model types
        for j = 1:numSteps
          mTypes{j} = cell([1 numOfStructures(j)]);
        end
        
        for j = 1:numSteps
            for k = 1:numOfStructures(j)
                % for lPar
                commonPart = {j,k,'lPar'};
                oneSturcturePar = [repmat(commonPart,numLPar,1),locsPar, cell(numLPar,1), repmat({""},numLPar,1)];
                allPars = [allPars; oneSturcturePar];
                
                % for mPar
                commonPart = {j,k,'mPar'};
                if contains(allModels(selMod{j}(k)),'.mat')
                    pars.par = [];
                    modelPar = cell(0,1);
                    numMPar = 0;
                    mTypes{j}{k} = 'image';
                else
                    modelfunc = str2func(allModels(selMod{j}(k))); % get model parameters
                    pars = modelfunc();
                    modelPar = fieldnames(pars.par);
                    numMPar = size(modelPar,1);
                    mTypes{j}{k} = pars.type;
                end
                oneSturcturePar = [repmat(commonPart,numMPar,1),modelPar, cell(numMPar,1), repmat({""},numMPar,1)];
                allPars = [allPars; oneSturcturePar];
                
            end
        end
        obj.setPar('allPars',allPars)
        obj.setPar('modelTypes',mTypes)
        assignSymbolOptions = cellstr([repmat('step', size(allPars,1),1) num2str(allPars.Step) repmat('.structure', size(allPars,1),1) num2str(allPars.Structure) repmat('.', size(allPars,1),1) char(allPars.ParType) repmat('.', size(allPars,1),1) char(allPars.Parameters)]);
        obj.guihandles.linkPar.String = assignSymbolOptions;
        obj.guihandles.linkPar.Enable = 'on';
        delete(gcf)
    end
end

function assign_callback(a,b,obj)
    allPars = obj.getPar('allPars');
    linkPar = obj.getSingleGuiParameter('linkPar');
    symbol = obj.getSingleGuiParameter('symbol');
    linkPar = linkPar.String(linkPar.Value);
    parParts = getParParts(linkPar);
    lLinkPar = allPars.Step==str2num(parParts(1)) & allPars.Structure==str2num(parParts(2)) & allPars.ParType==parParts(3) & allPars.Parameters==parParts(4);
    allPars.Symbol(lLinkPar) = symbol;
    obj.setPar('allPars', allPars);
end

function defineNewPar_callback(a,b,obj)
    % required info
    userDefinedPars = obj.getPar('userDefinedPars');
    expression = obj.getSingleGuiParameter('expression');
    symbol2 = obj.getSingleGuiParameter('symbol2');
    
    % save the symbo/expression pair
    lexisted = ismember(userDefinedPars.symbol,string(['userDefined.' symbol2]));
    if sum(lexisted)>0                  % if the same symbol exists
        userDefinedPars.expression(lexisted) = expression;
    else
        userDefinedPars = [userDefinedPars;{['userDefined.' symbol2] expression}];
    end
    obj.setPar('userDefinedPars', userDefinedPars);
    lPar = obj.getPar('lPar');
    mPar = obj.getPar('mPar');
    
    % update the options in lPar and mPar table
    for k=1:size(lPar,2)                % the kth step
        for l = 1:size(lPar{k},2)       % the lth structure
            lPar{k}{l}.Fix_fromFit = addcats(lPar{k}{l}.Fix_fromFit, userDefinedPars.symbol);
        end
    end
    for k=1:size(mPar,2)
        for l = 1:size(mPar{k},2)
            mPar{k}{l}.Fix_fromFit = addcats(mPar{k}{l}.Fix_fromFit, userDefinedPars.symbol);
        end
    end
    obj.setPar('lPar', lPar);
    obj.setPar('mPar', mPar);
end

function parParts = getParParts(x)
    parParts = string(strsplit(char(x),'.'));
    parParts(1) = strrep(parParts(1),'step','');
    parParts(2) = strrep(parParts(2),'structure','');
end
    
function ipTemplate = intrinsicParTemp()
    ipTemplate = [];
    ipTemplate.ch = [];
    ipTemplate.x = [];
    ipTemplate.y = [];
    ipTemplate.z = [];
    ipTemplate.xrot = [];
    ipTemplate.yrot = [];
    ipTemplate.zrot = [];
    ipTemplate.xscale = [];
    ipTemplate.yscale = [];
    ipTemplate.zscale = [];
    ipTemplate.weight = [];
end

function epTemplate = extrinsicParTemp()
    epTemplate = [];
    epTemplate.azimuthalShift = [];
    epTemplate.ringDistance = [];
    epTemplate.diameter = [];
end

function userDefinedPars = userDefinedParsTemp()
    userDefinedPars = table('Size',[0 2],'VariableTypes',{'string','string'},'VariableNames',{'symbol','expression'});
end