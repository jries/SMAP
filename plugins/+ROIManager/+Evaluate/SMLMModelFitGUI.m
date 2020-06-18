classdef SMLMModelFitGUI<interfaces.SEEvaluationProcessor
    properties
        fitter              % An SMLMModelFit object.
        numMod              % Number of component models.
        parsArgFieldnames   % 
        lFnParsArgEdit      %
        fnParsArgColWidth   %
        layerFieldnames     %
        lFnLayerEdit        %
        currentLoadedModel  %
        alignSettings       % Converter for alignment.
    end
    methods
        function obj=SMLMModelFitGUI(varargin)
            obj@interfaces.SEEvaluationProcessor(varargin{:});
            flagDirExist = exist('../SMLMModelFit','dir');
            if flagDirExist==0
                addpath(genpath('../ries-private'))
            else
                addpath(genpath('../SMLMModelFit'))
            end
            obj.propertiesToSave={'fitter', 'numMod', 'parsArgFieldnames', 'lFnParsArgEdit', 'fnParsArgColWidth', 'layerFieldnames', 'lFnLayerEdit', 'currentLoadedModel'};         
            addlistener(obj, 'mParsArgModified', @mParsArgModified_callback);
        end
        
        function setGuiParameters(obj,p)
            obj.fitter = p.fitter;
            initTabWhenLoading(obj);
            setGuiParameters@interfaces.SEEvaluationProcessor(obj,p);
        end
        
        function initTabWhenLoading(obj)
            fitter = obj.fitter;
            numOfModel = fitter.numOfModel;
            inp = obj.getAllParameters;
            obj.setPar('loading',true);
            for m = 1:numOfModel
                modelnumberStr = num2str(m);
                if ~isfield(inp,['modelname_' modelnumberStr])
                    % hack the action of pressing '+' tab
                    eventData.EventName = 'SelectionChanged';
                    eventData.OldValue = [];
                    eventData.NewValue = [];
                    eventData.NewValue.Title = '+';
                    obj.guihandles.tabgroup.SelectionChangedFcn{1}(obj.guihandles.tabgroup,eventData,obj)
                end
                % load info to the GUI
                obj.currentLoadedModel = modelnumberStr;
                notify(obj.guihandles.(['modelload_' modelnumberStr]),'Action');
            end
            obj.setPar('loading',false);
        end

        function out=run(obj, inp, varargin)
            % varargin is used only when called by the run_addguitotab.mlx
            p = inputParser;
            addParameter(p,'onlySetUp',false, @islogical);
            addParameter(p,'forceDisplay',false, @islogical);
            addParameter(p,'keepParsVal',false, @islogical);
            parse(p,varargin{:});
            results = p.Results;
            try
                out=runSMLMModelFitGUI(obj, inp, results.onlySetUp, results.forceDisplay, results.keepParsVal);
            catch
                warning(['Model fitter did not run through. Site ' num2str(obj.site.ID) ' encountered some issues.'])
                out.fitInfo = 'Fit failed.';
            end
        end
        
        function makeGui(obj,varargin)
            %% init settings
            % create the SMLMModelFit obj
                % check the dim of the data and update the corresponding
                % seeting in the fitter
                if isfield(obj.locData.loc,'znm')
                    dataDim = 3;
                else
                    dataDim = 2;
                end
            obj.fitter = SMLMModelFit('DataDim',dataDim);
            obj.fitter.linkedGUI = obj;
                
            obj.currentLoadedModel = [];
            
            % field names in the parsArg table
            fnParsArg={'name','value','fix','lb','ub','type','min','max','label'};
            lFnParsArgEdit = logical([0 1 1 1 1 0 1 1 1]);
            fnParsArgColWidth = {50,35,20,30,30,30,30,30,30};
            obj.parsArgFieldnames = fnParsArg;
            obj.lFnParsArgEdit = lFnParsArgEdit;
            obj.fnParsArgColWidth = fnParsArgColWidth;
            % field names in the layer table
            fnLayer ={'layer','value','fix','lb','ub','min','max'};
            fnLayerColWidth = {30,35,20,30,30,30,30};
            lFnLayerEdit = logical([0 1 1 1 1 1 1]);
            obj.layerFieldnames = fnLayer;
            obj.lFnLayerEdit = lFnLayerEdit;
            
            if nargin >2 && varargin{2}==true %gui of tabs
                 makeGui@interfaces.GuiModuleInterface(obj,varargin{1});
            else %real make GUI of main gui
                makeGui@interfaces.GuiModuleInterface(obj); %make the main GUI from the guidef definitions
                %Settings               
                optimizernames=obj.fitter.supportedSolver;
                obj.guihandles.optimizer.String=optimizernames;
                % Use fminsearchbnd by default.
                obj.guihandles.optimizer.Value = find(strcmp(optimizernames, 'fminsearchbnd'));
                Data={};
                oldh=obj.guihandles.optimizerpar;
                pos=oldh.Position;
                htable=uitable(oldh.Parent,'Data',Data,'Position',pos);
                htable.ColumnEditable = true;
                htable.ColumnName = {'Setting', 'Value'};
                htable.CellSelectionCallback = {@optimizerSetting_selectionCallback,obj};
                htable.CellEditCallback = {@optimizerSetting_editCallback,obj};
                htable.RowName = [];
                delete(oldh);
                obj.guihandles.optimizerpar=htable;
                notify(obj.guihandles.optimizer, 'Action');
                
                %% Layer settings
                oldh = obj.guihandles.layerSetting;
                pos=oldh.Position;
                hLayer=uitable(oldh.Parent,'Data',{},'Position',pos);
                hLayer.ColumnName = fnLayer;
                hLayer.ColumnEditable = true;
                hLayer.CellEditCallback = {@layerSetting_callback,obj};
                hLayer.RowName = [];
                hLayer.ColumnWidth = fnLayerColWidth;
                delete(oldh);
                obj.guihandles.layerSetting=hLayer;
                
                %% add first module and  + tab
                obj.guihandles.tabgroup=obj.guihandles.tab1.Parent;
                obj.addguitotab(1);
                obj.guihandles.addtab=uitab(obj.guihandles.tabgroup,'Title','+');
                obj.numMod = 1;                 % init of the model counts
                obj.guihandles.tabgroup.SelectionChangedFcn={@selectLayer_callback,obj};
                
                % Select the M1 tab by default since a user usually starts from loadin a model.
                obj.guihandles.tabgroup.SelectedTab = obj.guihandles.tab1;
                
                %% Converter tab 
                addconverttotab(obj);           % create the converter 
                oldh=obj.guihandles.anchorConvert;
                
                % create the converter input table
                htable = createConvertTable(oldh, obj);
                obj.guihandles.anchorConvert=htable;
                
                
            end
            obj.setPar('initiated',true)
        end
        function pard=guidef(obj)
            pard=SMLMModelFitGUIdef(obj);
        end
        
        function addguitotab(obj,number)
            run_addguitotab(obj,number);
        end
        
        function addconverttotab(obj)
            run_addconverttotab(obj);
        end
               
        function parId = loadParTable(obj, htable, fitter, modelnumber)
            % get parId and update the GUIParTable
            [parId,subParsArgTemp] = fitter.getAllParId(modelnumber);
            htable.Data = struct2Data(subParsArgTemp);
            htable.CellEditCallback = {@parSetting_callback,obj, modelnumber};
            htable.ColumnEditable = obj.lFnParsArgEdit;
            htable.ColumnWidth = obj.fnParsArgColWidth;
        end
        
        %%
        function mParsArgModified_callback(obj,b)
            modelbnumberStr = num2str(b.modelID);
            htable = obj.guihandles.(['partable_' modelbnumberStr ]);
            parId = obj.loadParTable(htable, obj.fitter, b.modelID);
            % here the target options in the converter table are defined.
            hConvert = obj.guihandles.anchorConvert;
            optionTarget = unique([hConvert.ColumnFormat{3} parId]);
            hConvert.ColumnFormat{3} = optionTarget;
            obj.guihandles.anchorConvert=hConvert;
        end
        
        function set.fitter(obj,value) 
            obj.fitter = value;
            if isfield(obj.P.par.mainGui.content.children.guiSites.children.Segment.processors,'SimulateSites')
                simulateSites = obj.P.par.mainGui.content.children.guiSites.children.Segment.processors.SimulateSites;
                if strcmp(simulateSites.guihandles.useFitter_button.Visible, 'off')
                    obj.P.par.mainGui.content.children.guiSites.children.Segment.processors.SimulateSites.initGui;
                end
            end
        end
    end
    events
        mParsArgModified
    end
end

%% Callbacks
function parSetting_callback(a,b,obj, modelnumber)
% Callback for editing parsArg of the model
% First get the fitter obj, and then reset the initial guess based on what
% are saved. This is to prevent any interference between different sites.
% Next detect the position of the modified value, check its counterpart in
% fitter.allParsArg, and overwrite it. Finally the initial guess is saved
% and the fitter of the obj replaced.

fitter = obj.fitter;
fitter.resetInit;
fn = obj.parsArgFieldnames;
indices = b.Indices;
data = a.Data;
name = strtrim(data{indices(1),strcmp(fn,'name')});
type = strtrim(data{indices(1),strcmp(fn,'type')});
[~,idx] = fitter.wherePar(['pars.m' num2str(modelnumber) '.' type '.' name]);
fitter.allParsArg.(fn{indices(2)})(idx) = b.NewData;
fitter.saveInit;
obj.fitter = fitter;
end

function selectLayer_callback(tabgroup,eventdata,obj) 
    % if + tab selected this makes a new model
    layertitle=(eventdata.NewValue.Title);
    if strcmp(layertitle,'+')
        numtabs=length(tabgroup.Children);
        obj.addguitotab(numtabs-2)
        obj.numMod = numtabs-2;   % save the number of model
        s=1:length(tabgroup.Children);
        % shift the order of table
        s(end-2)=s(end);
        s(end)=s(end)-2;
        s(end-1)=s(end);
        s(end)=s(end)+1;
        tabgroup.Children=tabgroup.Children(s);
        tabgroup.SelectedTab=tabgroup.Children(end-2); 
    end
end

function layerSetting_callback(a,b,obj) 
    fitter = obj.fitter;
    fn=obj.layerFieldnames;
    indices = b.Indices;
    layer = a.Data{indices(1),1};
    layer = strrep(layer,'layer','9');
    [~,idx] = fitter.wherePar(['pars.m' layer '.offset.weight']);
    fitter.allParsArg.(fn{indices(2)})(idx) = b.NewData;
    obj.fitter = fitter;
end

function optimizerSetting_selectionCallback(a,b,obj)
    try
    selectedRow = b.Indices(1);
    obj.setPar(['selectedRowOptimizer_' obj.name], selectedRow)
    catch
    end
end

function optimizerSetting_editCallback(a,b,obj)
    obj.setPar(['lOptimizerSettingEdited_' obj.name], true)
end