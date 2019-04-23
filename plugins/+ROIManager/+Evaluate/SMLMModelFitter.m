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
            pard=guidef(obj);
        end
    end
end

function pard = guidef(obj)
    pard.t1.object = struct('Style','text','String', 'Number of structures:');
    pard.t1.position = [1 1];
    pard.t1.Width = 2;
    
    pard.numOfStructure.object = struct('Style','popupmenu','String', 1:3, 'Value', 1);
    pard.numOfStructure.position = [1 3];
    pard.numOfStructure.Width = 0.5;
    
    pard.t2.object = struct('Style','text','String', 'Number of steps:');
    pard.t2.position = [2 1];
    pard.t2.Width = 2;
    
    pard.numOfStep.object = struct('Style','popupmenu','String', 1:3, 'Value', 1, 'Callback', {{@setStep_callback,obj}});
    pard.numOfStep.position = [2 3];
    pard.numOfStep.Width = 0.5;
    
    pard.currentStep.object = struct('Style','popupmenu','String', {'1'}, 'Value', 1);
    pard.currentStep.position = [2 3.5];
    pard.currentStep.Width = 0.5;
    
    pard.setIPar.object = struct('Style','pushbutton','String', 'Intrinsic Par', 'Callback', {{@setIPar_callback,obj}});
    pard.setIPar.position = [3 1];
    pard.setIPar.Width = 1;
    
    pard.setEPar.object = struct('Style','pushbutton','String', 'Extrinsic Par', 'Callback', {{@setEPar_callback,obj}});
    pard.setEPar.position = [3 2];
    pard.setEPar.Width = 1;
    
    pard.setModel.object = struct('Style','pushbutton','String', 'Select models', 'Callback', {{@selectModels_callback,obj}});
    pard.setModel.position = [3 3];
    pard.setModel.Width = 1;
    
    pard.displayFolder.object = struct('Style','edit','String', '---', 'Enable','off');
    pard.displayFolder.position = [4 1];
    pard.displayFolder.Width = 3;
    
    pard.modelFolder.object = struct('Style','pushbutton','String', 'Model', 'Callback', {{@load_callback,obj}});
    pard.modelFolder.position = [4 4];
    pard.modelFolder.Width = 1;
    
    pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi'};
    pard.plugininfo.type='ROI_Evaluate';
end

function setIPar_callback(obj,a,b)
    guiPar = b.getGuiParameters;
    numOfJ = guiPar.numOfStep.Value;
    numOfK = guiPar.numOfStructure.Value;
    k = str2double(guiPar.numOfStructure.selection);
    % k is the number of structures
    
    % GUI part
    f = figure('Name', 'Intrinsic_parameters');
    tabgp = uitabgroup(f,'Position',[.01 .01 0.98 0.98]);
    tab=cell(1:k);
    for l = 1:k
        tab{l} = uitab(tabgp,'Title',num2str(l));
    end
    
    % data part
    if isempty(b.getPar('iPar'))
        % initiate parameter table
        pars = fieldnames(intrinsicParTemp());
        used = num2cell(false(size(pars)));
        range = num2cell(zeros(size(pars)));
        oneData = [pars used range range used range];
        data=cell(numOfK,numOfJ);
        for l = 1:numOfK
            for m = 1:numOfJ
                data{l,m}=oneData;
            end
        end
        b.setPar('iPar', data);
    end
    currentK = 1; % this can be not just 1 later depends on the sets of parameters
    currentJ = guiPar.currentStep.Value;
    iPar = b.getPar('iPar');
	oneIPar = iPar{currentK, currentJ};
    uit = uitable(tab{currentK},'Data',oneIPar, 'Position', [20 20 400 204]);
    uit.ColumnName = {'Parameter','Fix','Lower bound','Upper bound','Specify int', 'Int'};
    uit.ColumnEditable = [false true true true true true];
    uit.CellEditCallback = @updatePars;
    
    function updatePars(src,event)
        iPar{currentK, currentJ} = uit.Data;
        b.setPar('iPar',iPar)
    end
end

function setEPar_callback(obj,a,b)
    guiPar = b.getGuiParameters;
    k = str2double(guiPar.numOfStructure.selection);
    f = figure('Name', 'Extrinsic_parameters');
    tabgp = uitabgroup(f,'Position',[.01 .01 0.98 0.98]);
    tab=cell(1:k);
    for l = 1:k
        tab{l} = uitab(tabgp,'Title',num2str(l));
    end
    currentK = 1; % this can be not just 1 later depends on the sets of parameters
    if isempty(b.getPar('ePar'))
        % initiate parameter table
        pars = fieldnames(extrinsicParTemp());
        used = num2cell(false(size(pars)));
        range = num2cell(zeros(size(pars)));
        oneData = [pars used range range used range];
        data={};
        data{1}=oneData;
        b.setPar('ePar', data);
    end
    ePar = b.getPar('ePar');
	oneEPar = ePar{currentK};
    uit = uitable(tab{currentK},'Data',oneEPar, 'Position', [20 20 400 204]);
    uit.ColumnName = {'Parameter','Fix','Lower bound','Upper bound','Specify int', 'Int'};
    uit.ColumnEditable = [false true true true true true];
    uit.CellEditCallback = @updatePars;
    
    function updatePars(src,event)
        ePar{currentK} = uit.Data;
        b.setPar('ePar',ePar)
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
    allModels = ls([f '\*.m']);
    obj.setGuiParameters(struct('displayFolder',f));
    obj.setPar('allModels',allModels);
end

function setStep_callback(a,b,obj)
    guiPars = obj.getGuiParameters;
    guiPars.currentStep.String = num2str((1:a.Value)');
    obj.setGuiParameters(guiPars)
end

function selectModels_callback(a,b,obj)
    f = figure('Name', 'Select Models');
    selMod = obj.getPar('selectedModel');
    c = {};
    if isempty(selMod)
        selMod = {};
    else
    end
    numOfK = obj.getSingleGuiParameter('numOfStructure').Value;
    numOfJ = obj.getSingleGuiParameter('numOfStep').Value;
    for l = 1:numOfJ
        for m = 1:numOfK
            if isempty(selMod)|~isequal(size(selMod),[numOfJ numOfK])
                selMod{l,m} = 1;
            end
            c{l,m} = uicontrol(f, 'Style','popupmenu','String',obj.getPar('allModels'),'Value',selMod{l,m});
            c{l,m}.Position = [20+(l-1)*210 350+(m-1)*20 200 20];
            c{l,m}.Callback = {@updatePars,{obj,l,m}};
        end
    end
    
    function updatePars(a,b,obj)
        selMod{obj{2},obj{3}} = b.Source.Value;
        obj{1}.setPar('selectedModel',selMod);
    end
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