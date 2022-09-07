classdef LocMoFitGUI<interfaces.SEEvaluationProcessor
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
        sourceModel         % The source function model of the image model.
        compiledMode        %
    end
    properties (Dependent)
        lPreview            % which model is for previewing
    end
    methods
        function obj=LocMoFitGUI(varargin)
            obj@interfaces.SEEvaluationProcessor(varargin{:});
            obj.propertiesToSave={'fitter', 'numMod', 'parsArgFieldnames', 'lFnParsArgEdit', 'fnParsArgColWidth', 'layerFieldnames', 'lFnLayerEdit', 'currentLoadedModel','sourceModel'};
            addlistener(obj, 'mParsArgModified', @mParsArgModified_callback);
        end
               
        function setGuiParameters(obj,p,setchildren,setmenu)
            obj.fitter = p.fitter;
            obj.fitter.updateVersion;
            initTabWhenLoading(obj);
            % For upgrade from SMLMModelFit to LocMoFit
            if ~isempty(p.anchorConvert.Data)
                col_source = p.anchorConvert.Data(:,1);
                col_source = replace(col_source, 'SMLMModelFitGUI', 'LocMoFitGUI');
                p.anchorConvert.Data(:,1) = col_source;
            end
            
            % Check the model type are consistent between the GUI and obj.
            m = 1;
            while isfield(p, ['modelType_' num2str(m)])
                ID = ['modelType_' num2str(m)];
                p.(ID).String = obj.guihandles.(ID).String;
                p.(ID).Value = obj.guihandles.(ID).Value;
                if iscell(p.(ID).String)
                    p.(ID).selection = p.(ID).String{p.(ID).Value};
                else
                    p.(ID).selection = p.(ID).String;
                end
                m = m+1;
            end
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
            out = run_LocMoFitGUI(obj, inp, varargin{:});
        end
        
        function makeGui(obj,varargin)
            %% init settings
            % create the SMLMModelFit obj
            % check the dim of the data and update the corresponding
            % setting in the fitter
            if isfield(obj.locData.loc,'znm')
                dataDim = 3;
            else
                dataDim = 2;
            end
            obj.fitter = SMLMModelFit('DataDim',dataDim);
            
            obj.currentLoadedModel = [];

            isdeployed = true;
            if isdeployed
                obj.compiledMode = 'on';
                model2Load = modelList;
                hModels = [];
                for k = 1:length(model2Load)
                    tempModel = modelList(model2Load{k});
                    hModels.(model2Load{k}) = tempModel;
                end
                obj.setPar('modelOptions', hModels)
            else
                obj.compiledMode = 'off';
            end
            
            %% ParArg table related things
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
                if ismac
                    obj.guiPar.FieldHeight=20;
                    obj.guiPar.fontsize=12;
                    obj.guiPar.Xrim=20;
                    obj.guiPar.Vrim=0;
                    obj.guiPar.Vpos=1.5;
                    obj.guiPar.Xpos=0.9;
%                     tabh=[50 160];
                end
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
                htable.ColumnName = {'Parameters', 'Value'};
                htable.ColumnWidth = {100, 60};
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
                
                %% add first module and + tab
                obj.guihandles.tabgroup=obj.guihandles.tab1.Parent;
                obj.addguitotab(1);
                obj.guihandles.rmtab=uitab(obj.guihandles.tabgroup,'Title','[x]');
                obj.guihandles.addtab=uitab(obj.guihandles.tabgroup,'Title','+');
                obj.numMod = 1;                 % init of the model counts
                obj.guihandles.tabgroup.SelectionChangedFcn={@selectLayer_callback,obj};
                
                % Select the M1 tab by default since a user usually starts from loading a model.
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
            if ismac
                dy=1;
            else
                dy=0;
            end

            pard.tab.tab1='Settings';
            
            pard.t_optimizer.object=struct('Style','text','String','Optimizer:');
            pard.t_optimizer.position=[1.5+dy,1];
            pard.t_optimizer.Width=1;
            pard.t_optimizer.tab='tab1';
            
            pard.optimizer.object=struct('Style','popupmenu','String',{{'text','t2'}}, 'Callback',{{@optimizer_callback,obj}});
            pard.optimizer.position=[1.5+dy,2];
            pard.optimizer.Width=1.5;
            pard.optimizer.tab='tab1';
            pard.optimizer.TooltipString = 'Select the optimizer for the fitting here.';
            
            pard.loadfitting.object=struct('Style','pushbutton','String','load', 'Callback',{{@obj.load_callback}});
            pard.loadfitting.position=[1+dy,3.8];
            pard.loadfitting.Width=0.6;
            pard.loadfitting.Height=0.7;
            pard.loadfitting.tab='tab1';
            pard.loadfitting.TooltipString='Load previously saved settings.';
            
            pard.savefitting.object=struct('Style','pushbutton','String','save', 'Callback',{{@save_callback,obj}});
            pard.savefitting.position=[1+dy,4.4];
            pard.savefitting.Width=0.6;
            pard.savefitting.Height=0.7;
            pard.savefitting.tab='tab1';
            pard.savefitting.TooltipString='Save the current settings.';
            
            pard.setting_advanced.object=struct('Style','pushbutton','String','Advanced','Callback',{{@setting_advanced_callback,obj}});
            pard.setting_advanced.position=[1.8+dy,3.8];
            pard.setting_advanced.Width=1.2;
            pard.setting_advanced.Height=0.7;
            pard.setting_advanced.tab='tab1';
            pard.setting_advanced.TooltipString='Advanced settings.';
            
            pard.noFit.object=struct('Style','checkbox','String','Review only','Value',0);
            pard.noFit.position=[3.5+dy,3.5];
            pard.noFit.Width=1.3;
            pard.noFit.Height=1;
            pard.noFit.tab='tab1';
            pard.noFit.TooltipString='If checked, the fit will not be excecuted. This is useful for viewing the old fit result.';
                       
            pard.useAlignment.object=struct('Style','checkbox','String','Register','Value',0);
            pard.useAlignment.position=[4.5+dy,3.5];
            pard.useAlignment.Width=1.5;
            pard.useAlignment.Height=1;
            pard.useAlignment.tab='tab1';
            pard.useAlignment.TooltipString = 'Regester the current sites to a common coordinate system. If needed, click on [...] to define extra transformations before the registration.';
            
            pard.setting_alignment.object=struct('Style','pushbutton','String','...','Callback',{{@setting_alignment_callback,obj}});
            pard.setting_alignment.position=[4.3+dy,4.7];
            pard.setting_alignment.Width=0.25;
            pard.setting_alignment.Height=0.5;
            pard.setting_alignment.tab='tab1';
            pard.setting_alignment.TooltipString = 'Settings for extra transformations.';
                                  
            pard.optimizerpar.object=struct('Style','text','String','');
            pard.optimizerpar.position=[4.5+dy,1];
            pard.optimizerpar.Width=2.5;
            pard.optimizerpar.Height=3;
            pard.optimizerpar.tab='tab1';
            pard.optimizerpar.TooltipString = 'Define settings for the optimizer of choice here.';
            
            pard.addRowOptimizer.object=struct('Style','pushbutton','String','+', 'Callback',{{@addRowOptimizer_callback,obj}});
            pard.addRowOptimizer.position=[5+dy,3];
            pard.addRowOptimizer.Width=0.2;
            pard.addRowOptimizer.Height=0.4;
            pard.addRowOptimizer.tab='tab1';
            pard.addRowOptimizer.TooltipString = 'Add a new row.';
            
            pard.rmRowOptimizer.object=struct('Style','pushbutton','String','-', 'Callback',{{@rmRowOptimizer_callback,obj}});
            pard.rmRowOptimizer.position=[5+dy,3.2];
            pard.rmRowOptimizer.Width=0.2;
            pard.rmRowOptimizer.Height=0.4;
            pard.rmRowOptimizer.tab='tab1';
            pard.rmRowOptimizer.TooltipString = 'Remove the selected row.';
            
            pard.t_layerSetting.object=struct('Style','text','String','Layer background:');
            pard.t_layerSetting.position=[6+dy,1];
            pard.t_layerSetting.Width=1.5;
            pard.t_layerSetting.tab='tab1';
            
            pard.layerSetting.object=struct('Style','text','String','');
            pard.layerSetting.position=[10+dy,1];
            pard.layerSetting.Width=3.9;
            pard.layerSetting.Height=4;
            pard.layerSetting.tab='tab1';

            pard.helpButton.object=struct('Style','pushbutton','String','?', 'Callback',{{@helpLocMoFitGUI,obj}});
            pard.helpButton.position=[0+dy,4.7];
            pard.helpButton.Width=0.3;
            pard.helpButton.Height=0.6;
            pard.helpButton.tab='none';
            pard.helpButton.TooltipString = 'Help page for the current tab.';
        end
        
        function addguitotab(obj,number)
            %% RUN_ADDGUITOTAB Adding a model (number) to the tab group.
            %
            tag=['M' num2str(number)];
            obj.guihandles.(['tab' num2str(number)])=uitab(obj.guihandles.tabgroup,'Title',tag,'Tag',tag);
            pardModelTab = guidefmodel(obj,number);
            pardModelTab = obj.setVisibility_compiledMode(pardModelTab);
            guidefhere=addnumbertofield(pardModelTab,number);
            guiParold=obj.guiPar;handleold=obj.handle;
            obj.guiPar.Vrim=0;
            if ismac
                obj.guiPar.FieldHeight=20;
                obj.guiPar.fontsize=12;
                obj.guiPar.Xrim=25;
                obj.guiPar.Vrim=0;
                obj.guiPar.Xpos=0.7;
                obj.guiPar.Vpos=1.5;
                tabh=[50 160];
            else
                tabh=[20 120];
            end
            obj.handle=obj.guihandles.(['tab' num2str(number)]);
            fitter = obj.fitter;
            obj.makeGui(guidefhere,1);
            obj.fitter = fitter;
            obj.fitter.linkedGUI = obj;
            obj.handle=handleold;
            obj.guiPar=guiParold;
            %initialize parameter table. Here change the looks!
            hpar=obj.guihandles.(['tabpar_' num2str(number)]);
            hsettings=obj.guihandles.(['tabsettings_' num2str(number)]);
            
            % parameter tab
            
            % flatten the guihandles
            obj.guihandles.(['partable_' num2str(number)])=uitable(hpar,'Data',{});
            obj.guihandles.(['partable_' num2str(number)]).ColumnName = obj.parsArgFieldnames;
            obj.guihandles.(['partable_' num2str(number)]).Position(3:4)=obj.guihandles.(['partable_' num2str(number)]).Position(3:4)-tabh;
            obj.guihandles.(['partable_' num2str(number)]).ColumnEditable = obj.lFnParsArgEdit;
            obj.guihandles.(['partable_' num2str(number)]).RowName = [];
            
            % button for saving the model
            obj.guihandles.(['savePar_' num2str(number)])=uicontrol(hpar,'Style','pushbutton','String','Export');
            obj.guihandles.(['savePar_' num2str(number)]).Position = [20 0 40 20];
            obj.guihandles.(['savePar_' num2str(number)]).Callback = {@savePar_callback, obj, number};
            obj.guihandles.(['savePar_' num2str(number)]).Tooltip = 'Export the settings as a text file.';

            % button for loading the model
            obj.guihandles.(['loadPar_' num2str(number)])=uicontrol(hpar,'Style','pushbutton','String','Import');
            obj.guihandles.(['loadPar_' num2str(number)]).Position = [60 0 40 20];
            obj.guihandles.(['loadPar_' num2str(number)]).Callback = {@loadPar_callback, obj, number};
            obj.guihandles.(['loadPar_' num2str(number)]).Tooltip = 'Import previously exported settings.';
            
            % set settings as default
            obj.guihandles.(['setDefPar_' num2str(number)])=uicontrol(hpar,'Style','pushbutton','String','set def.', 'Visible', 'off');
            obj.guihandles.(['setDefPar_' num2str(number)]).Position = [100 0 35 20];
            obj.guihandles.(['setDefPar_' num2str(number)]).Callback = {@setDefPar_callback, obj};
            obj.guihandles.(['setDefPar_' num2str(number)]).Tooltip = 'Set the current settings as default.';
            
            % write all settings
            obj.guihandles.(['writePar_' num2str(number)])=uicontrol(hpar,'Style','pushbutton','String','write', 'Visible', 'off');
            obj.guihandles.(['writePar_' num2str(number)]).Position = [135 0 35 20];
            obj.guihandles.(['writePar_' num2str(number)]).Callback = {@writePar_callback, obj};
            obj.guihandles.(['writePar_' num2str(number)]).Tooltip = 'Write the settings for the current site.';
            
            % read all settings
            obj.guihandles.(['readPar_' num2str(number)])=uicontrol(hpar,'Style','pushbutton','String','read', 'Visible', 'off');
            obj.guihandles.(['readPar_' num2str(number)]).Position = [170 0 35 20];
            obj.guihandles.(['readPar_' num2str(number)]).Callback = {@readPar_callback, obj};
            obj.guihandles.(['readPar_' num2str(number)]).Tooltip = 'Read the previously written settings to the GUI.';
            
            % clear all settings
            obj.guihandles.(['clearPar_' num2str(number)])=uicontrol(hpar,'Style','pushbutton','String','reset', 'Visible', 'off');
            obj.guihandles.(['clearPar_' num2str(number)]).Position = [205 0 35 20];
            obj.guihandles.(['clearPar_' num2str(number)]).Callback = {@clearPar_callback, obj};
            obj.guihandles.(['clearPar_' num2str(number)]).Tooltip = 'Clear the written settings for the current site.';
            
            % checkbox for previewing the model
            obj.guihandles.(['forceDisplay_' num2str(number)])=uicontrol(hpar,'Style','checkbox','String','Preview', 'Value',0);
            obj.guihandles.(['forceDisplay_' num2str(number)]).Position = [240 0 120 20];
            obj.guihandles.(['forceDisplay_' num2str(number)]).Tooltip = 'When this option is checked, the fit will be skipped and the model is visualized based on the starting paramters.';
%             obj.guihandles.(['forceDisplay_' num2str(number)]).Callback = {@forceDisplay_callback, obj, number};
            
            % internal settings tab
            obj.guihandles.(['settingstable_' num2str(number)])=uitable(hsettings,'Data',{});
            obj.guihandles.(['settingstable_' num2str(number)]).ColumnName = {'setting','value'};
            obj.guihandles.(['settingstable_' num2str(number)]).Position(3:4)=obj.guihandles.(['settingstable_' num2str(number)]).Position(3:4)-tabh;
            obj.guihandles.(['settingstable_' num2str(number)]).ColumnEditable = true;
            obj.guihandles.(['settingstable_' num2str(number)]).RowName = [];
            obj.guihandles.(['settingstable_' num2str(number)]).CellSelectionCallback = {@settingTable_cellSelectionCallback, obj, num2str(number)};
            
%             % Buttons for adding removing a row of setting
%             obj.guihandles.(['addInternalSetting_' num2str(number)])=uicontrol(hsettings,'Style','pushbutton','String','+');
%             obj.guihandles.(['addInternalSetting_' num2str(number)]).Position = [260 0 20 20];
%             obj.guihandles.(['addInternalSetting_' num2str(number)]).Callback = {@addInternalSetting_callback, obj, num2str(number)};
%             obj.guihandles.(['rmInternalSetting_' num2str(number)])=uicontrol(hsettings,'Style','pushbutton','String','-');
%             obj.guihandles.(['rmInternalSetting_' num2str(number)]).Position = [280 0 20 20];
%             obj.guihandles.(['rmInternalSetting_' num2str(number)]).Callback = {@rmInternalSetting_callback, obj, num2str(number)};
        end

        function rmguifromtab(obj,number)
            %% rmguifromtab removes the last model from the GUI.
            fn = fieldnames(obj.guihandles);
            indRm = find(endsWith(fn,['_' num2str(number)]));
            for k = length(indRm):-1:1
                delete(obj.guihandles.(fn{indRm(k)}));
                obj.guihandles = rmfield(obj.guihandles, fn{indRm(k)});
            end
        end
        
        function addconverttotab(obj)
            tag='Convert';
            obj.guihandles.converter=uitab(obj.guihandles.tabgroup,'Title',tag,'Tag',tag);
            guiParold=obj.guiPar;handleold=obj.handle;
            obj.guiPar.Vrim=0;
            obj.handle=obj.guihandles.converter;
            if ismac
                obj.guiPar.Vrim = 65;
                obj.guiPar.Vpos = 0;
                obj.guiPar.Xrim = 20;
                obj.guiPar.Xpos = 0.9;
            else
                obj.guiPar.Vrim = 45;
                obj.guiPar.FieldHeight = 22;
            end
            guidefhere=guidefconvert(obj);
            obj.makeGui(guidefhere,1);
            obj.handle=handleold;
            obj.guiPar=guiParold;
        end
        
        function parId = loadParTable(obj, htable, fitter, modelnumber)
            % get parId and update the GUIParTable
            [parId,subParsArgTemp] = fitter.getAllParId(modelnumber, 'form', 'long');
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
        
        %% advanced settings
        function setAdvanceSetting(obj, par, value)
            obj.advanceSetting.(par) = value;
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
        %% Updating the layer settings for the GUI
        function updateLayer(obj)
            % Update the layers in use.
            % Get the fitter.
            fitter = obj.fitter;
            hLayer=obj.guihandles.layerSetting; % Hande of the layer table.
            % Get the current settings for the layers.
            layerParsArg = fitter.subParsArg([]);
            % Convert cellstr into chr array so that uitalble can handle.
            fn = obj.layerFieldnames;
            layerStr = num2str(layerParsArg.model);
            layerParsArgDisp.layer = char(strcat('layer', layerStr(:,2)));
            for k = 2:length(fn)
                if iscell(layerParsArg.(fn{k}))
                    layerParsArgDisp.(fn{k}) = char(layerParsArg.(fn{k}));
                else
                    layerParsArgDisp.(fn{k}) = layerParsArg.(fn{k});
                end
            end
            hLayer.Data = struct2Data(layerParsArgDisp);
            hLayer.ColumnEditable = obj.lFnLayerEdit;
            % save the table back to the GUI
            obj.guihandles.layerSetting = hLayer;
            layerParId = {};
            for k = length(layerParsArg.model):-1:1
                layerParId(k) = fitter.getAllParId(layerParsArg.model(k));
            end
            %% update the convert table
            hConvert = obj.guihandles.anchorConvert;
            colFormat = hConvert.ColumnFormat{3};
            % remove old layer parID
            lLayerParID = contains(colFormat,'.offset.weight');
            if any(lLayerParID)
                colFormat(lLayerParID) = [];
            end
            % add new layer parID
            optionTarget = unique([colFormat layerParId{:}]);
            hConvert.ColumnFormat{3} = optionTarget;
            obj.guihandles.anchorConvert=hConvert;
        end

        % Configure the visibility in response to certain change here
        function parModelTab = setVisibility_compiledMode(obj, parModelTab)
            visible_whenOn = {'modelname_CM','modelload_CM'};
            visible_whenOff = {'modelname','modelload'};
            switch obj.compiledMode
                case 'on'
                    visible = visible_whenOn;
                    invisible = visible_whenOff;
                case 'off'
                    visible = visible_whenOff;
                    invisible = visible_whenOn;
            end
            for k = 1:length(visible)
                parModelTab.(visible{k}).Visible = 'on';
            end
            for k = 1:length(invisible)
                parModelTab.(invisible{k}).Visible = 'off';
            end
        end

        % callback for loading
        function load_callback(obj,a,b,path2setting)
            % get the path of the saved info
            if obj.getPar('initiated')
                if ~exist('path2setting', 'var')
                    [file,path2setting] = uigetfile({'*_SMLMModelFit.mat;*_LocMoFit.mat','LocMoFit settings'}, 'Select an SMLMModelFit file...', '');
                else
                    [path2setting,file,ext] = fileparts(path2setting);
                    file = [file ext];
                    path2setting = [path2setting filesep];
                end
                temp = load([path2setting file]);
                obj.setGuiParameters(temp.par);
                obj.fitter = temp.par.fitter;
                obj.fitter.linkedGUI = obj;
                obj.fitter.loadListener;
                for k = 1:obj.fitter.numOfModel
                    obj.fitter.model{k}.loadListener;
                end
            end
        end

        function updateGUI_fromLocMoFitObj(obj)
            fitter = obj.fitter;
            for m = 1:fitter.numOfModel
                htable = obj.guihandles.(['partable_' num2str(m)]);
                obj.loadParTable(htable, fitter, m);
                obj.updateAdvanceTab(fitter, m);
                obj.setPar('loading',true)
                initmodel(obj, num2str(m),'skipAddModel',true, 'compiledMode', false);
                obj.setPar('loading',false)
            end

            for m = 1:fitter.numOfModel
                modelStr = num2str(m);
                obj.guihandles.(['setDefPar_' modelStr]).Visible = fitter.getAdvanceSetting('siteSpecificMode');
                obj.guihandles.(['writePar_' modelStr]).Visible = fitter.getAdvanceSetting('siteSpecificMode');
                obj.guihandles.(['clearPar_' modelStr]).Visible = fitter.getAdvanceSetting('siteSpecificMode');
                obj.guihandles.(['readPar_' modelStr]).Visible = fitter.getAdvanceSetting('siteSpecificMode');
            end
        end

        function updateGUI_convert_fromLocMoFitObj(obj)
            % Uupdate GUI: update the convert table based on the LocMoFit
            % Obj.
            % Known issue: targets start with 'usr_' cannot be processed.
            fitter = obj.fitter;
            allSourceFitter = fitter.converterSource;
            loaded_beforeMe = obj.find_LocMoFit_beforeMe;
            fn = fieldnames(loaded_beforeMe);
            fitterGUI_order = {};
            for k = 2:length(allSourceFitter)
                for m = 1:length(fn)
                    l = isequal(allSourceFitter{k}, loaded_beforeMe.(fn{m}));
                    if l
                        fitterGUI_order{k} = fn{m};
                    end
                end
            end
            numOfRules = length(fitter.converterRules.target);
            data = cell(numOfRules,4);
            data(:,1:3) = [fitterGUI_order(fitter.converterRules.target_Id)' fitter.converterRules.rule_raw' fitter.converterRules.target'];
            obj.guihandles.anchorConvert.Data = data;
        end

        function updateAdvanceTab(obj,fitter, modelnumberStr)
            m = modelnumberStr;
            htable = obj.guihandles.(['settingstable_' num2str(m)]);
            internalPar_list = fitter.getModelInternalSettingList(m);
            data = {};
            for k = length(internalPar_list):-1:1
                currentPar = internalPar_list{k};
                val = fitter.getModelInternalSetting(m, currentPar);
                data(k,:) = {currentPar val};
            end
            htable.Data = data;
        end

        function fitter_stack = find_LocMoFit_beforeMe(obj)
            se = obj.locData.SE;
            eval_ = se.processors.eval;
            loadedPlugins = fieldnames(eval_.children);
            loadedLocMoFitGUI = loadedPlugins(startsWith(loadedPlugins, 'LocMoFit'));
            currentName = obj.name;
            ind_current = find(ismember(loadedLocMoFitGUI, currentName));
            if ind_current > 1
                loadedLocMoFitGUI_beforeMe = loadedLocMoFitGUI(1:ind_current-1);
                fn = fieldnames(se.processors.eval.children);
                for k = 1:length(loadedLocMoFitGUI_beforeMe)
                    fitter_stack.(loadedLocMoFitGUI_beforeMe{k}) = se.processors.eval.children.(loadedLocMoFitGUI_beforeMe{k}).fitter;
                end
            else
                fitter_stack = [];
            end
            
        end
        function value = get.lPreview(obj)
            inp = obj.getGuiParameters;
            fn = fieldnames(inp);
            ind = startsWith(fieldnames(inp),'forceDisplay_');
            fn_subset = fn(ind);
            forceDisplay = [];
            for k = length(fn_subset):-1:1
                forceDisplay(k) = inp.(fn_subset{k});
            end
            value = forceDisplay;
        end
        function resetPar(obj)
            obj.fitter.resetInit
            obj.fitter.setAllModelInternalSetting(obj.getPar('allModIntSetting_default'));
            obj.guihandles.anchorConvert.Data = obj.getPar('data_convert');
        %     obj.updateGUI_convert_fromLocMoFitObj
            obj.updateGUI_fromLocMoFitObj
        end
        function readPar(obj)
            allModIntSetting = obj.site.evaluation.(obj.name).written.allModIntSetting;
            obj.fitter.setAllModelInternalSetting(allModIntSetting)      %IC220907 -Moved 3 lines up  
            obj.fitter.allParsArg = obj.site.evaluation.(obj.name).written.allParsArg;
            obj.guihandles.anchorConvert.Data = obj.site.evaluation.(obj.name).written.data_convert;
            
        
        %     obj.updateGUI_convert_fromLocMoFitObj
            obj.updateGUI_fromLocMoFitObj
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

siteSpecificMode = fitter.getAdvanceSetting('siteSpecificMode');
switch siteSpecificMode
    case 'off'
        fitter.saveInit;
        fitter.lockInit;
    case 'on'
%         Do nothing
end

obj.fitter = fitter;
end

function selectLayer_callback(tabgroup,eventdata,obj)
% if + tab selected this makes a new model
layertitle=(eventdata.NewValue.Title);
switch layertitle
    case '+'
        numberOfExtraTab = 3;
        numtabs=length(tabgroup.Children);
        obj.addguitotab(numtabs-numberOfExtraTab)
        obj.numMod = numtabs-numberOfExtraTab;   % save the number of model
        s=1:length(tabgroup.Children);
        % shift the order of table
        s = [s(1:end-numberOfExtraTab-1) s(end) s(end-numberOfExtraTab:end-1)];
        tabgroup.Children=tabgroup.Children(s);
        tabgroup.SelectedTab=tabgroup.Children(end-numberOfExtraTab);
    case '[x]'
        allTabTitle = {tabgroup.Children.Title};
        idxRmButton = find(strcmp(layertitle, allTabTitle));
        mod2rm = allTabTitle{idxRmButton - 1};
        rmDone = 0;
        if strcmp(mod2rm, 'M1')
            obj.setPar('status', 'this model cannot be removed (see warning)')
            warning(['There should be at least one model so ' mod2rm ' cannot be removed.'])
        else
            answer = questdlg(['Remove ' mod2rm '?'],...
                'Remove the last model', ...
                'Yes', 'No',...
                'No');
            switch answer
                case 'Yes'
                    modelID = str2num(mod2rm(2:end));
                    obj.rmguifromtab(modelID);
                    delete(obj.guihandles.(['tab' mod2rm(2:end)]));
                    obj.guihandles = rmfield(obj.guihandles, ['tab' mod2rm(2:end)]);
                    if obj.fitter.numOfModel == modelID
                        obj.fitter.rmLastModel;
                    end
                    obj.numMod = obj.numMod-1;
                    obj.setPar('status',[mod2rm ' has been successfully removed'])
                    rmDone = 1;
                otherwise

            end
            
        end
        tabgroup.SelectedTab = tabgroup.Children(idxRmButton-rmDone-1);
    otherwise
end
end

function layerSetting_callback(a,b,obj)
fitter = obj.fitter;
fn=obj.layerFieldnames;
indices = b.Indices;
layer = a.Data{indices(1),1};
layer = strrep(layer,'layer','9');
offsetForm = obj.fitter.getAdvanceSetting(['m' num2str(layer) '_background']);
[~,idx] = fitter.wherePar(['pars.m' layer '.offset.' offsetForm]);
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

%% Format conversion
% From uitableData to parsArg struct
function parsArg = data2ParsArg(uitableData,fn,model)
parsArgOri = cell2struct(uitableData,fn(2:end),2);
for k = 2:length(fn)
    if ischar([parsArgOri.(fn{k})])
        tempField = strtrim({parsArgOri.(fn{k})});
    else
        tempField = [parsArgOri.(fn{k})];
    end
    parsArg.(fn{k}) = tempField';
end
parsArg.model = repelem(model, length(parsArg.(fn{2})))';
parsArg = orderfields(parsArg, fn);
end

%% guidef related


% Call back for saving SMLMModelFit settings

function save_callback(a,b,obj)
[file,path] = uiputfile('*_LocMoFit.mat', 'Save as', '');
obj.fitter.resetInit;
gm = obj.getPar('mainGui');
ge = gm.children.guiSites.children.eval;
pe = ge.getGuiParameters(true);
par = pe.children.(obj.name);
var2save = {'par'};
save([path file], [var2save{:}],'-v7.3');
end

%         Create a new row in the optimizer table

function addRowOptimizer_callback(a,b,obj)
htable = obj.guihandles.optimizerpar;
htable.Data = [htable.Data; {[],[]}];
obj.guihandles.optimizerpar=htable;
end

%         Delete a row in the optimizer table

function rmRowOptimizer_callback(a,b,obj)
htable = obj.guihandles.optimizerpar;
data = htable.Data;
selectedRow = obj.getPar(['selectedRowOptimizer_'  obj.name]);
fitter = obj.fitter;
try
    data(selectedRow,:) = [];
catch
end
htable.Data = data;
obj.guihandles.optimizerpar=htable;
end

%         Call back for selecting the optimizer

% function forceDisplay_callback(a,b,obj,modelID)
%     currentName = obj.name;
%     GUINumber = strrep(currentName, 'LocMoFitGUI', '');
%     if isempty(GUINumber)
%         GUINumber = 1;
%     else
%         GUINumber = str2num(GUINumber(2:end));
%     end
%     LocMoFitGUI_preview = obj.locData.getPar('LocMoFitGUI_preview');
%     LocMoFitGUI_preview(GUINumber, modelID) = a.Value;
%     obj.locData.setPar('LocMoFitGUI_preview', LocMoFitGUI_preview);
% end

function optimizer_callback(a,b,obj)
% !!! later here should update the options for the optimizer settings
% based on the optimizer selected.
selectedOptimizer = b.Source.String{b.Source.Value};
hSetting = obj.guihandles.optimizerpar;
hSetting.ColumnFormat = {getSolverOption(selectedOptimizer),[]};
obj.setPar(['lOptimizerSettingEdited_' obj.name], true);
end

function setting_advanced_callback(a,b,obj)
fig = figure('Name','Advance setting');
advanceSetting = obj.fitter.advanceSetting;
fn = fieldnames(advanceSetting);
nrow = length(fn);
fig.Position(3:4) = [320 nrow*22+60];
guihandles = [];
for k = 1:nrow
    option = advanceSetting.(fn{k}).option;
    value = advanceSetting.(fn{k}).value;
    name = advanceSetting.(fn{k}).name;
    t = ['uit_' num2str(k)];
    ui = ['ui_' num2str(k)];
    lineHeight = 0.73;
    if islogical(value)
        guihandles.(t) = uicontrol(fig,'Style','checkbox','String',name,'Value',value);
        guihandles.(t).Position = [1.2 (k-1)*lineHeight 1.5 0.9];
    else
        if length(option)>1
            idxVal = find(strcmp(value,option));
            guihandles.(t) = uicontrol(fig,'Style','text','String',name);
            guihandles.(t).Position = [1.2 (k-1)*lineHeight 1.5 0.9];
            guihandles.(ui) = uicontrol(fig,'Style','popupmenu','Value',idxVal,'String',option);
            guihandles.(ui).Position = [2.7 (k-1)*lineHeight 1 0.9];
        else
            guihandles.(t) = uicontrol(fig,'Style','text','String',name);
            guihandles.(t).Position = [1.2 (k-1)*lineHeight 1.5 0.9];
            guihandles.(ui) = uicontrol(fig,'Style','edit','String',num2str(value));
            guihandles.(ui).Max = 2;
            guihandles.(ui).Position = [2.7 (k-1)*lineHeight 1 0.9];
        end
    end
    
    if islogical(value)
        obj.guihandles.([tag 't_' fn{k}]) = guihandles.(t);
    end
    
    obj.guihandles.(fn{k}) = guihandles.(ui);
end

guihandles.save = uicontrol(fig,'Style','pushbutton','String','Save','Callback',{@saveAdvanceSettings_callback, fig, obj});
guihandles.save.Position = [3.2 k*lineHeight+0.3 0.5 0.9];

guiStyle(guihandles, fieldnames(guihandles),'FieldHeight',20);
end

function saveAdvanceSettings_callback(a,b,fig,obj)
fitter = obj.fitter;
advanceSetting = fitter.advanceSetting;
fn = fieldnames(advanceSetting);
for k = 1:length(fn)
    oneCtrl = obj.guihandles.(fn{k});
    switch oneCtrl.Style
        case 'popupmenu'
            fitter.setAdvanceSetting(fn{k},oneCtrl.String{oneCtrl.Value});
        case 'edit'
            fitter.setAdvanceSetting(fn{k},str2num(oneCtrl.String));
        case 'checkbox'
            fitter.setAdvanceSetting(fn{k},oneCtrl.String);
    end
end
try
    obj.updateGUI_fromLocMoFitObj
    obj.updateLayer
catch
    warning('Layer(s) is not updated.')
end
close(fig)
end

%         Call back for alignment settings

function setting_alignment_callback(a,b,obj)
fig = figure(514);
fig.Name = 'Transformation settings';
fig.Position(3:4) = [400 360];
guihandles.uit = uitable(fig);
guihandles.uit = createConvertTable(guihandles.uit, fig);
guihandles.uit.CellEditCallback = {@alignSettingEdit_callback,3,obj};

% Update the convert tab.
parId = {'post_x', 'post_y', 'post_z', 'post_scale', 'post_zrot'};
optionTarget = unique([guihandles.uit.ColumnFormat{2} parId]);
guihandles.uit.ColumnFormat{2} = optionTarget;
guihandles.uit.Data = obj.alignSettings;
guihandles.uit.ColumnWidth = {100 100 100};



guihandles.button_add = uicontrol('Style','pushbutton','String','+','Callback',{@addNewRuleAlign_callback,fig,obj});
guihandles.button_add.Position = [1 10 0.3 0.8];

guihandles.button_rm = uicontrol('Style','pushbutton','String','-','Callback',{@rmRuleAlign_callback,fig,obj});
guihandles.button_rm.Position = [1.3 10 0.3 0.8];
guihandles.uit.Position = [1 0 3.7 10];
guiStyle(guihandles, fieldnames(guihandles));
end

function addNewRuleAlign_callback(a,b,fig,obj)
htable = findobj(fig, 'Type', 'uitable');
htable.Data = [htable.Data; {[],[],[]}];
end

function rmRuleAlign_callback(a,b,fig,obj)
htable = findobj(fig, 'Type', 'uitable');
data = htable.Data;
temp = findobj(fig,'Type','uicontrol','-and','String','selectedRowConvert');
try
    selectedRow = temp.Value;
    data(selectedRow,:) = [];
    htable.Data = data;
    obj.alignSettings = htable.Data;
catch
    display('Please select a row first.')
end

end

function alignSettingEdit_callback(a,b,k,obj)
obj.alignSettings = a.Data;
end

%% addguitotab related
function out=addnumbertofield(in,number)
% A helper function. obj.guihandles needs to be flat. Add number to field name.
fn=fieldnames(in);
for k=1:length(fn)
    if strcmp(fn{k},'tab')
        out.tab=addnumbertofield(in.(fn{k}),number);
    else
        out.([fn{k} '_' num2str(number)])=in.(fn{k});
    end
end
end

% GUI for the model tab
function pard=guidefmodel(obj,number)
if ismac
    dy=1.5;
else
    dy=0;
end
pard.tab.tabmodel='Model';

pard.modelname.object=struct('Style','edit','String','');
pard.modelname.position=[3+dy,1];
pard.modelname.Width=1.8;
pard.modelname.tab=['tabmodel_' num2str(number)];
pard.modelname.Visible = 'on';

pard.modelload.object=struct('Style','pushbutton','String','load model','Callback',{{@loadmodel_callback,obj}});
pard.modelload.position=[3+dy,3];
pard.modelload.Width=1;
pard.modelload.tab=['tabmodel_' num2str(number)];
pard.modelload.Visible = 'on';

fn = {};
modelOptions = obj.getPar('modelOptions');
if ~isempty(modelOptions)
    fn = fieldnames(modelOptions);
end
fn = ['Select a model...'; fn; '[from a file...]'];
pard.modelname_CM.object=struct('Style','popupmenu','String', {fn},'value', 1);
pard.modelname_CM.position=[3+dy,1];
pard.modelname_CM.Width=1.8;
pard.modelname_CM.tab=['tabmodel_' num2str(number)];
pard.modelname_CM.Visible = 'off';
pard.modelname_CM.Tooltip = 'The geometric model to be loaded.';

pard.modelInfo.object=struct('Style','pushbutton','String','i','Callback',{{@helpLocMoFitGUI,obj,true}});
pard.modelInfo.position=[3+dy,2.8];
pard.modelInfo.Width=0.2;
pard.modelInfo.tab=['tabmodel_' num2str(number)];
pard.modelInfo.Tooltip = 'Information of the selected model.';

pard.modelload_CM.object=struct('Style','pushbutton','String','load model','Callback',{{@loadmodel_callback,obj}});
pard.modelload_CM.position=[3+dy,3];
pard.modelload_CM.Width=1;
pard.modelload_CM.tab=['tabmodel_' num2str(number)];
pard.modelload_CM.Visible = 'off';
pard.modelload_CM.Tooltip = 'Load the selected model to LocMoFit.';

pard.modelType.object=struct('Style','popupmenu','String',{'Type'},'value', 1,'Callback',{{@modType_callback,obj,number}});
pard.modelType.position=[3+dy,4];
pard.modelType.Width=1;
pard.modelType.tab=['tabmodel_' num2str(number)];
pard.modelType.Tooltip='Convert the model to the type you specify.';
pard.modelType.Enable = 'off';
pard.modelType.Tooltip = 'Model forms.';

pard.layert.object=struct('Style','text','String','Layer');
pard.layert.position=[4+dy,1];
pard.layert.Width=2;
pard.layert.tab=['tabmodel_' num2str(number)];
pard.layert.Tooltip = 'The layer that the model is fit to.';

pard.layer.object=struct('Style','popupmenu','String',{{'1','2','3','4','5','6'}}, 'Callback', {{@layer_callback,obj,number}});
pard.layer.position=[4+dy,3.5];
pard.layer.Width=1;
pard.layer.tab=['tabmodel_' num2str(number)];
pard.layer.Tooltip = "This corresponds to the layer defined in the 'Render' tab.";
pard.layer.Enable = 'off';

pard.pixelsizefitt.object=struct('Style','text','String','Pixel size');
pard.pixelsizefitt.position=[5+dy,1];
pard.pixelsizefitt.Width=2;
pard.pixelsizefitt.tab=['tabmodel_' num2str(number)];
pard.pixelsizefitt.Tooltip = 'The bin size for binning the model.';

pard.pixelsizefit.object=struct('Style','edit','String','10','Callback',{{@modelSettingsEdited_callback,obj}});
pard.pixelsizefit.position=[5+dy,3.5];
pard.pixelsizefit.Width=1;
pard.pixelsizefit.tab=['tabmodel_' num2str(number)];
pard.pixelsizefit.Tooltip = 'The smaller the more precise.';
pard.pixelsizefit.Enable = 'off';

pard.sigma_fitt.object=struct('Style','text','String','Sigma of Gaussian filtering');
pard.sigma_fitt.position=[6+dy,1];
pard.sigma_fitt.Width=2;
pard.sigma_fitt.tab=['tabmodel_' num2str(number)];
pard.sigma_fitt.Tooltip = 'The Gaussian filtering makes the model smoother.';

pard.sigma_fit.object=struct('Style','edit','String','15','Callback',{{@modelSettingsEdited_callback,obj}});
pard.sigma_fit.position=[6+dy,3.5];
pard.sigma_fit.Width=1;
pard.sigma_fit.tab=['tabmodel_' num2str(number)];
pard.sigma_fit.Enable = 'off';

pard.sigmaOrFactor.object=struct('Style','text','String','');
pard.sigmaOrFactor.position=[6+dy,4.5];
pard.sigmaOrFactor.Width=1;
pard.sigmaOrFactor.Height=2;
pard.sigmaOrFactor.tab=['tabmodel_' num2str(number)];
pard.sigmaOrFactor.Tooltip = 'Select the mode of sigma to use.';

pard.useSigma.object=struct('Style','radiobutton', 'Value', 1,'Callback',{{@modelSettingsEdited_callback,obj}});
pard.useSigma.position=[1+dy,1];
pard.useSigma.Width=1;
pard.useSigma.tab=['tabmodel_' num2str(number)];
pard.useSigma.Visible = 'off';

pard.sigmaFactor_fitt.object=struct('Style','text','String','Factor of Gaussian filtering');
pard.sigmaFactor_fitt.position=[7+dy,1];
pard.sigmaFactor_fitt.Width=2;
pard.sigmaFactor_fitt.tab=['tabmodel_' num2str(number)];
pard.sigmaFactor_fitt.Tooltip = 'The Gaussian filtering makes the model smoother.';

pard.useSigmaFactor.object=struct('Style','radiobutton','Value', 0,'Callback',{{@modelSettingsEdited_callback,obj}});
pard.useSigmaFactor.position=[2+dy,1];
pard.useSigmaFactor.Width=1;
pard.useSigmaFactor.tab=['tabmodel_' num2str(number)];
pard.useSigmaFactor.Visible = 'off';

pard.sigmaFactor_fit.object=struct('Style','edit','String','1','Callback',{{@modelSettingsEdited_callback,obj}});
pard.sigmaFactor_fit.position=[7+dy,3.5];
pard.sigmaFactor_fit.Width=1;
pard.sigmaFactor_fit.tab=['tabmodel_' num2str(number)];
pard.sigmaFactor_fit.Enable = 'off';

% the sub tabs for the model are created here
pard.tab.tabpar='Parameters';
pard.tab.tabsettings='Advance';
end

% Actions when loading a model
function loadmodel_callback(a,b,obj)
% Executed when loading a model.
% Load model based on the model number (ID).
currentLoadedModel = obj.currentLoadedModel;
if isempty(currentLoadedModel)
    modelnumber=(a.Parent.Parent.Parent.Title(2:end)); %hack to get the right tab
    if ~isempty(obj.guihandles.(['modelname_CM_' modelnumber]).String)
        modelOptions = obj.guihandles.(['modelname_CM_' modelnumber]).String;
        modelSelected = modelOptions{obj.guihandles.(['modelname_CM_' modelnumber]).Value};
    else
        modelSelected = '';
    end
    if strcmp(obj.compiledMode,'off')||startsWith(modelSelected, '[')
        filter = {'*.m;*.mlx;*_img.mat;*.png','Supported formats (**.m,*.mlx,*_img.mat,*.png)';...
            '*.m;*.mlx','Code files (*.m,*.mlx)'; ...
            '*_img.mat','Matrix files (*_img.mat)'; ...
            '*.png','Image files (*.png)'; ...
            };
        if startsWith(modelSelected, '[')
            if strcmp(modelSelected, '[from a file...]')
                fnold = '';
            else
                fnold=modelSelected(2:end-1);
            end
        else
            fnold = obj.guihandles.(['modelname_' modelnumber]).String;
        end
        if isempty(fnold)
            [~,fnold] = LocMoFit.getModelList;
        end
        [f,p]=uigetfile(filter,'Select a geometric model',fnold);
        if ~f %no model selected: return
            obj.setPar('status','no model loaded')
            disp('Cancelled by the user. No model loaded.')
            return
        end
        if startsWith(modelSelected, '[')
            toAdd = ['[' p f ']'];
            idxToAdd = find(strcmp(toAdd, obj.guihandles.(['modelname_CM_' modelnumber]).String));
            if ~isempty(idxToAdd)
                obj.guihandles.(['modelname_CM_' modelnumber]).Value = idxToAdd;
            else
                obj.guihandles.(['modelname_CM_' modelnumber]).String{end+1}=['[' p f ']'];
                obj.guihandles.(['modelname_CM_' modelnumber]).Value = length(obj.guihandles.(['modelname_CM_' modelnumber]).String);
            end
        else
            obj.guihandles.(['modelname_' modelnumber]).String=[p f];
        end
    end
else
    modelnumber = currentLoadedModel;
end
% If the model is loaded from a saved file, then skip addModel and load the
% model obj from the file, otherwise create a new obj.
if obj.getPar('loading')
    init_status = initmodel(obj, modelnumber,'skipAddModel',true, 'compiledMode', obj.compiledMode);
else
    init_status = initmodel(obj, modelnumber, 'compiledMode', obj.compiledMode);
end


% Update the layers in use.
if init_status
    obj.updateLayer;
    obj.setPar('status', 'model successfully loaded')
else
    obj.setPar('status', 'no model is loaded: please select a model first')
end
end

%    Initiate the modelwhich
function status = initmodel(obj, modelnumber,varargin)
% Initiate the model.
% Initiate the function.
p = inputParser;
addParameter(p, 'skipAddModel',false)
addParameter(p, 'compiledMode','off')
parse(p,varargin{:});
results = p.Results;
fitter = obj.fitter;
% Get model's path.
modPath = obj.guihandles.(['modelname_' modelnumber]).String;
if ~strcmp(modPath,'Loaded')&&~results.skipAddModel
    
    %% Decide which subclass of the SMLMModel to use based on the input model.
    switch results.compiledMode
        case 'off'
            lFromFile = true;
            fileSource = modPath;
        case 'on'
            modelOptions = obj.guihandles.(['modelname_CM_' modelnumber]).String;
            idxModelSeleted = obj.guihandles.(['modelname_CM_' modelnumber]).Value;
            modelSeleted = modelOptions{idxModelSeleted};
            if startsWith(modelSeleted,'[')
                lFromFile = true;
                fileSource = modelSeleted(2:end-1);
                modPath = fileSource;
            elseif endsWith(modelSeleted,'...')
                disp('Please select a model first.')
                status = 0;
                return
            else
                lFromFile = false;
            end
    end
    
    if lFromFile
        [filePath, modelFun, ext] = fileparts(fileSource);

        % Add the root of the model script as a path

        %     modelFun = str2func(modelFun);
        %     tempModelObj = modelFun();
        if ismember(ext,{'.png', '.bmp', '.tif', '.tiff'}) || (strcmp(ext,{'.mat'}) && endsWith(modelFun,'_img'))
            if strcmp(ext, '.mat')
                img = load(modPath,'img');
                img = img.img;
            else
                img = imread(modPath);
            end
            geoModeltemp = imageModel(img);
        else
            addpath(filePath);
            geoModeltemp = functionModel(modPath);
        end
    else
        oneGuiHandle = obj.guihandles.(['modelname_CM_' modelnumber]);
        ind = oneGuiHandle.Value;
        modelName = oneGuiHandle.String{ind};
        modelOptions = obj.getPar('modelOptions');
        geoModeltemp = functionModel(copy(modelOptions.(modelName)));
    end

    if isempty(obj.locData.loc)
        % When there is not localizations, assuming the user is loading the
        % geometric model for simulations.
        fitter.dataDim = geoModeltemp.dimension;
    end

    % If the position of the model has been occupied, replace the occupying model.
    if ~isempty(fitter.model)&&~length(fitter.model)<str2double(modelnumber)
        fitter.changeModel(geoModeltemp, str2double(modelnumber));
    else
        fitter.addModel(geoModeltemp);
    end

end
% Set up the GUI of the model tab.
enableSettings(obj,modelnumber,fitter);
modelnumberStr = modelnumber;
modelnumber = str2double(modelnumber);
% Load model's info to the GUI.
obj.guihandles.(['layer_' modelnumberStr]).Value = fitter.model{modelnumber}.layer;
if ismember(fitter.model{modelnumber}.modelType,{'discrete','continuous'})
    obj.guihandles.(['sigma_fit_' modelnumberStr]).String = fitter.model{modelnumber}.sigma;
    obj.guihandles.(['sigmaFactor_fit_' modelnumberStr]).String = num2str(fitter.model{modelnumber}.sigmaFactor);
end
obj.guihandles.(['pixelsizefit_' modelnumberStr]).String = fitter.model{modelnumber}.pixelSize;

%             _for the parameter table_
% Update the model parameter table when the model has been loaded.
htable=obj.guihandles.(['partable_' modelnumberStr]); %handle of table
parId = obj.loadParTable(htable, fitter, modelnumber);
obj.fitter = fitter;

%             _for the internal setting table (in the advanced tab)_
% Assign the callback function and default options.
hSettingTable=obj.guihandles.(['settingstable_' modelnumberStr]); %handle of table
hSettingTable.CellEditCallback = {@settingTable_cellEditCallback,obj, modelnumber};

modelInternalSetting = fitter.getModelInternalSettingList(modelnumber);
if isempty(modelInternalSetting)
    hSettingTable.ColumnFormat{1} = [];
    hSettingTable.Data = {};
else
    % !!!200528: I was asked to input a row vector instead of column. Keep this in
    % mind if things happend again.
    hSettingTable.Data = {};
    hSettingTable.Data(:,1) = cellstr(modelInternalSetting)';
    for k = length(modelInternalSetting):-1:1
        hSettingTable.Data{k,2} = fitter.getModelInternalSetting(modelnumber,modelInternalSetting{k});
    end
    hSettingTable.ColumnEditable = [false true];
    hSettingTable.ColumnWidth = {150 50};
%     hSettingTable.
end

% Update the convert tab.
hConvert = obj.guihandles.anchorConvert;
optionTarget = unique([hConvert.ColumnFormat{3} parId]);
hConvert.ColumnFormat{3} = optionTarget;
obj.guihandles.anchorConvert=hConvert;

% Update the model type options.
    switch obj.fitter.model{modelnumber}.modelType
        case 'image'
            if isempty(obj.sourceModel)
                modelOption ={'image'};
            else
                modelOption ={obj.sourceModel{modelnumber}.modelObj.modelTypeOption{:} 'image'};
            end
        case 'locsImg'
            modelOption ={obj.fitter.model{modelnumber}.modelObj.modelTypeOption{:}};
        otherwise
            modelOption ={obj.fitter.model{modelnumber}.modelObj.modelTypeOption{:} 'image'};
    end
    indType = find(strcmp(modelOption,obj.fitter.model{modelnumber}.modelType));
    obj.guihandles.(['modelType_' num2str(modelnumber)]).String = modelOption;
    obj.guihandles.(['modelType_' num2str(modelnumber)]).Value = indType;

% Update the model options
if results.skipAddModel
    % update the path to the model class
    SMLMModelObj = obj.fitter.model{str2num(modelnumberStr)};
    if ~isempty(SMLMModelObj.modelObj)
        path_model = SMLMModelObj.sourcePath;
        modelClass = class(SMLMModelObj.modelObj);
    else
        path_model = 'Image saved in the object';
        modelClass = char();
    end
    
    modelOptions = fieldnames(obj.getPar('modelOptions'));
    [l, ind] = ismember(modelClass,modelOptions);
    if l
        obj.guihandles.(['modelname_CM_' modelnumberStr]).String = modelOptions;
        obj.guihandles.(['modelname_CM_' modelnumberStr]).Value = ind;
    else
        obj.guihandles.(['modelname_CM_' modelnumberStr]).String = [modelOptions;['[' path_model ']']];
        obj.guihandles.(['modelname_CM_' modelnumberStr]).Value = length(modelOptions)+1;
    end
end
status = 1;
end

% Enable settings
function enableSettings(obj,modelnumber,fitter)
% Enable settings according to the model type.
modelType = fitter.model{str2double(modelnumber)}.modelType;
if isequal(modelType, 'image')
    % for image models
    obj.guihandles.(['layer_' modelnumber]).Enable = 'on';
    obj.guihandles.(['pixelsizefit_' modelnumber]).Enable = 'on';
    obj.guihandles.(['sigma_fit_' modelnumber]).Enable = 'off';
    obj.guihandles.(['sigmaFactor_fit_' modelnumber]).Enable = 'off';
elseif isequal(modelType, 'continuous')||isequal(modelType, 'background')
    % for continuous models
    obj.guihandles.(['layer_' modelnumber]).Enable = 'on';
    obj.guihandles.(['pixelsizefit_' modelnumber]).Enable = 'on';
    obj.guihandles.(['sigma_fit_' modelnumber]).Enable = 'on';
    obj.guihandles.(['sigmaFactor_fit_' modelnumber]).Enable = 'off';
else
    % for discrete/discretized models
    obj.guihandles.(['layer_' modelnumber]).Enable = 'on';
    obj.guihandles.(['pixelsizefit_' modelnumber]).Enable = 'on';
    obj.guihandles.(['pixelsizefit_' modelnumber]).String = 1;
    
    % for point models, user can choose to set sigma or sigma factor
    % here create the UI for the selection based on radiobutton
    obj.guihandles.(['useSigmaFactor_' modelnumber]).Visible = 'on';
    obj.guihandles.(['useSigma_' modelnumber]).Visible = 'on';
    hParent = obj.guihandles.(['sigmaOrFactor_' modelnumber]).Parent;
    obj.guihandles.(['sigmaOrFactor_' modelnumber]) = uibuttongroup('BorderType','none');
    obj.guihandles.(['sigmaOrFactor_' modelnumber]).SelectionChangedFcn = {@sigmaOrFactor_SeChangedFcn ,obj, modelnumber};
    obj.guihandles.(['sigmaOrFactor_' modelnumber]).Position = [0.89 0.33 0.1 0.24];
    obj.guihandles.(['sigmaOrFactor_' modelnumber]).Parent = hParent;
    obj.guihandles.(['useSigmaFactor_' modelnumber]).Parent = obj.guihandles.(['sigmaOrFactor_' modelnumber]);
    obj.guihandles.(['useSigmaFactor_' modelnumber]).Position = [0 32 15 15]; %! the 32 here is used as an flag for the callback 'sigmaOrFactor_SeChangedFcn'
    obj.guihandles.(['useSigma_' modelnumber]).Parent = obj.guihandles.(['sigmaOrFactor_' modelnumber]);
    obj.guihandles.(['useSigma_' modelnumber]).Position = [0 7 15 15];

    if obj.guihandles.(['useSigma_' modelnumber]).Value
        obj.guihandles.(['sigma_fit_' modelnumber]).Enable = 'off';
        obj.guihandles.(['sigmaFactor_fit_' modelnumber]).Enable = 'on';
    else
        obj.guihandles.(['sigma_fit_' modelnumber]).Enable = 'on';
        obj.guihandles.(['sigmaFactor_fit_' modelnumber]).Enable = 'off';
    end
end
obj.guihandles.(['modelType_' modelnumber]).Enable = 'on';

obj.updateLayer;
end

% Callback functions
% Callback for checking a new layer
function layer_callback(a, b, obj, modelNumber)
fitter = obj.fitter;
fitter.resetInit;
fitter.model{modelNumber}.layer = a.Value;
fitter.updateLayer;
fitter.saveInit;
obj.fitter = fitter;
obj.updateLayer;
obj.setPar(['logicalParTableEdited_' obj.name], true); % flag for editing the model settings
end

function modType_callback(a, b, obj, modelNumber)
% Callback for the conversion to an image model.
% Usually only obj.modelType has to be changed unless it is changing
% between image and funciton types.
fitter = obj.fitter;
switch a.String{a.Value}
    case 'image'
        if strcmp(fitter.model{modelNumber}.modelType, 'image')
            warning('The model is already the image type.')
        else
            selection = questdlg('Convert model?','Convert the model to the image type? (Before you convert the model, please make sure you set the right parameters...)','Yes','No','No');
            if strcmp(selection,'Yes')
                prompt = {'Pixel size'};
                dlgtitle = 'Set the pixel size';
                dims = [1 35];
                definput = {'5'};
                pxSize = inputdlg(prompt,dlgtitle,dims,definput);
                pxSize = str2num(pxSize{1});
                % export the image based on the current parameter values and use it
                % as an image model
                
                % save the source model of the image model for later.
                obj.sourceModel{modelNumber} = fitter.model{modelNumber};
                modelImg = fitter.getImage(modelNumber, 'pixelSize', pxSize);
                imgModel = imageModel(modelImg, 'pixelSize', pxSize);
                fitter.changeModel(imgModel,modelNumber);
            end
        end
    otherwise
        if strcmp(fitter.model{modelNumber}.modelType, 'image')
            % do the following if the previous model form is 'image'
            funModel = obj.sourceModel{modelNumber};
            fitter.changeModel(funModel,modelNumber);
            obj.sourceModel{modelNumber} = [];
        else
            fitter.model{modelNumber}.modelType = a.String{a.Value};
        end
end
enableSettings(obj,num2str(modelNumber),fitter);
initmodel(obj, num2str(modelNumber),'skipAddModel',true);
end

% Callback for the model settings edited
function modelSettingsEdited_callback(a,b,obj)
obj.setPar(['logicalParTableEdited_' obj.name], true);
end

% Callback for settingTable
function settingTable_cellEditCallback(a,b,obj, modelnumber)
indices = b.Indices;
data = a.Data;
obj.fitter.linkedGUI = obj;
obj.fitter.loadListener;
for k = 1:obj.fitter.numOfModel
    obj.fitter.model{k}.loadListener;
end
if indices(2)==1
    setting = strtrim(data{indices(1),1});
    value = obj.fitter.getModelInternalSetting(modelnumber,setting);
    data{indices(1),2} = value;
else
    setting = strtrim(data{indices(1),1});
    value = data{indices(1),2};
    obj.fitter.setModelInternalSetting(modelnumber,setting,value);
end
obj.guihandles.(['settingstable_' num2str(modelnumber)]).Data = data;
obj.setPar(['selectedRow_internalSetting_' obj.name '_' num2str(modelnumber)],indices(1));
end

% Callback for settingTable
function settingTable_cellSelectionCallback(a,b,obj, modelnumber)
indices = b.Indices;
if ~isempty(indices)
    obj.setPar(['selectedRow_internalSetting_' obj.name '_' num2str(modelnumber)],indices(1));
end
end

% Callback for internal settings
function addInternalSetting_callback(a,b,obj, modelNumber)
htable = obj.guihandles.(['settingstable_' modelNumber]);
htable.Data = [htable.Data; {[],[]}];
obj.guihandles.(['settingstable_' modelNumber])=htable;
end

function rmInternalSetting_callback(a,b,obj, modelNumber)
htable = obj.guihandles.(['settingstable_' modelNumber]);
data = htable.Data;
selectedRow= obj.getPar(['selectedRow_internalSetting_' obj.name '_' modelNumber]);
try
    data(selectedRow,:) = [];
catch
end
htable.Data = data;
obj.guihandles.(['settingstable_' modelNumber]) = htable;
end

% Action when switching model type
function sigmaOrFactor_SeChangedFcn(a,b,obj, modelNumber)
if b.Source.SelectedObject.Position(2)==32
    obj.fitter.model{str2double(modelNumber)}.fixSigma = true;
    obj.guihandles.(['sigma_fit_' modelNumber]).Enable = 'on';
    obj.guihandles.(['sigmaFactor_fit_' modelNumber]).Enable = 'off';
else
    obj.fitter.model{str2double(modelNumber)}.fixSigma = false;
    obj.guihandles.(['sigma_fit_' modelNumber]).Enable = 'off';
    obj.guihandles.(['sigmaFactor_fit_' modelNumber]).Enable = 'on';
end
end

function savePar_callback(a,b,obj,modID)
% Saving the current parameter table
filter = {'*_fitPar.csv'};
[file, path] = uiputfile(filter);
parTable = obj.guihandles.(['partable_' num2str(modID)]).Data;
writecell(parTable, [path file],'Delimiter','tab')
end

function loadPar_callback(a,b,obj,modID)
% Loading and matching a parameter table
filter = {'*_fitPar.csv'};
[file, path] = uigetfile(filter);
importTable = readcell([path file],'Delimiter','\t','LineEnding','\r\n');
dim = size(importTable);
if dim(2) == 8
    importTable = [importTable repmat({char('')},dim(1),1)];
    importTable = cell2table(importTable);
    importTable.Properties.VariableNames = obj.guihandles.(['partable_' num2str(modID)]).ColumnName;
end
refTable = obj.guihandles.(['partable_' num2str(modID)]);
[keyStr,locR,locT] = matchTable(refTable,importTable,{'type','name'},'.');
fig = figure('Name','Import parameter settings');
fig.Position(3:4) = [300 420];
uit = uitable(fig);
uit.Position = [20 40 300-40 420-80];
uit.ColumnName = {'import', 'to'};
ID_tar = strcat(importTable.type,'.',importTable.name);
ID_ref = strcat(refTable.Data(:,6),'.',refTable.Data(:,1));
col2 = cell(size(ID_tar));
col2(locT) = ID_ref(locR);
lEmpty = cellfun(@isempty, col2);
numNotMatched = sum(lEmpty);
if numNotMatched>0
    col2(lEmpty) = cellstr(repmat('none', [numNotMatched 1]));
end
uit.Data = [ID_tar col2];
uit.ColumnFormat = {[],[ID_ref' {'none'}]};
uit.ColumnEditable = [false true];
uit.ColumnWidth = {100, 100};
button = uicontrol(fig,'Style','pushbutton', 'String','Import', 'Callback', {@importPar_callback,obj,modID,uit, importTable, ID_ref});
end

function importPar_callback(a,b,obj,modID,uit, importTable, ID_ref)
refTable = obj.guihandles.(['partable_' num2str(modID)]);
[~,lInput,lOri] = intersect(uit.Data(:,2),ID_ref);
imported_found = table2cell(importTable(lInput,[2:5 7:9]));
imported_found(:,2) = cellfun(@logical,imported_found(:,2), 'UniformOutput', false);
refTable.Data(lOri,[2:5 7:9]) = imported_found;
obj.fitter.setParArgBatch(refTable.Data, 'modelID', modID);
close(uit.Parent)
end

function writePar_callback(a,b,obj)
    % Write the current settings as the default for the current site.
    allModIntSetting = obj.fitter.getAllModelInternalSetting;
    obj.site.evaluation.(obj.name).written = [];
    obj.site.evaluation.(obj.name).written.allParsArg = obj.fitter.allParsArg;
    obj.site.evaluation.(obj.name).written.allModIntSetting = allModIntSetting;
    obj.site.evaluation.(obj.name).written.data_convert = obj.guihandles.anchorConvert.Data;

    obj.setPar('status', 'The current settings are successfully wirtten.');
end

function readPar_callback(a,b,obj)
    obj.readPar;
end

function clearPar_callback(a,b,obj)
    % Clear the written settings.
    answer = questdlg('If you select yes, your previous written settings will be lost and the current parameter settings will be overwritten by the default',...
        'Cear the written settings?', ...
        'Yes', 'No',...
        'No');
    switch answer
        case 'Yes'
            if isempty(obj.site)
                obj.setPar('status', 'Nothing to clear: no site is selected. Plesae select one site first.');
            elseif isfield(obj.site.evaluation, obj.name)&&isfield(obj.site.evaluation.(obj.name), 'written')
                obj.site.evaluation.(obj.name) = rmfield(obj.site.evaluation.(obj.name), 'written');
                obj.setPar('status', 'The current settings are successfully cleared.');
            else
                obj.setPar('status', 'Nothing to clear: no setting is wirtten.');
            end
            if ~isempty(obj.getPar('allModIntSetting_default'))
                resetPar(obj)
            else
            end
        otherwise
            obj.setPar('status', 'Nothing is cleared: user cancelled.');
    end
end

function setDefPar_callback(a,b,obj)
    obj.fitter.saveInit;
    obj.fitter.lockInit;
    obj.setPar('allModIntSetting_default', obj.fitter.getAllModelInternalSetting)
    obj.setPar('data_convert', obj.guihandles.anchorConvert.Data)
end

%

function [keyStr,locR,locT] = matchTable(ref,tar,keys, delimiter)
switch class(ref)
    case 'matlab.ui.control.Table'
        refTab = cell2table(ref.Data);
        refTab.Properties.VariableNames  = ref.ColumnName;
    case 'table'
        refTab = ref;
    otherwise
        error(['ref class "' class(ref) '" not supported. Please use uitable or table instead.'])
        return
end
for k = 1:length(keys)
    if k == 1
        keyStr_ref = refTab.(keys{k});
        keyStr_tar = tar.(keys{k});
    else
        keyStr_ref = strcat(keyStr_ref, delimiter, refTab.(keys{k}));
        keyStr_tar = strcat(keyStr_tar, delimiter, tar.(keys{k}));
        keyStr_ref = deblank(keyStr_ref);
        keyStr_tar = deblank(keyStr_tar);
    end
end
[keyStr,locR,locT] = intersect(keyStr_ref,keyStr_tar);
end

%% # addconverttotab

% ## Define the gui
function pard=guidefconvert(obj)
pard.anchorConvert.object=struct('Style','text','String','');
pard.anchorConvert.position=[9.5,1];
pard.anchorConvert.Width=4;
pard.anchorConvert.Height=9.5;

pard.addNewRule.object=struct('Style','pushbutton','String','+','Callback',{{@addNewRule_callback,obj,'anchorConvert'}});
pard.addNewRule.position=[10.1,1];
pard.addNewRule.Width=0.2;
pard.addNewRule.Height=0.4;

pard.rmRule.object=struct('Style','pushbutton','String','-','Callback',{{@rmRule_callback,obj,'anchorConvert'}});
pard.rmRule.position=[10.1,1.2];
pard.rmRule.Width=0.2;
pard.rmRule.Height=0.4;

pard.evokeMatchPar.object=struct('Style','pushbutton','String','Match','Callback',{{@evokeParMatcher_callback,obj}});
pard.evokeMatchPar.position=[10.5,4.2];
pard.evokeMatchPar.Width=0.8;
pard.evokeMatchPar.Height=0.8;
end

% ## Match parameters
function hStack = evokeParMatcher_callback(a,b,obj)
obj.setPar('matchPar_fitGUIselected',[]);
fig = figure;
sourceSMLMFit_t = uicontrol(fig, 'Style', 'text','Position',[10 380 130 20], 'String', 'Source LocMoFitGUI','HorizontalAlignment','left');
sourceSMLMFit = uicontrol(fig, 'Style', 'popupmenu','Position',[140 380 100 20]);
sourceSMLMFit.String = obj.guihandles.anchorConvert.ColumnFormat{1};

sourceModel_t = uicontrol(fig, 'Style', 'text','Position',[10 350 130 20], 'String', 'Source model','HorizontalAlignment','left');
sourceModel = uicontrol(fig, 'Style', 'popupmenu','Position',[140 350 100 20]);
sourceModel.String = {'--'};

matchTable_t = uicontrol(fig, 'Style', 'text','Position',[10 320 350 20], 'String', 'Assign to the selected parameters in this LocMoFitGUI:','HorizontalAlignment','left');

sourceSMLMFit.Callback = {@SMLMFitSelect_callback,obj, sourceModel};

matchTable = uitable(fig,'Position',[10 40 250 270]);
allPossibleTargets = obj.guihandles.anchorConvert.ColumnFormat{3};
matchTable.Data = [num2cell(true(size(allPossibleTargets))); allPossibleTargets]';
matchTable.ColumnEditable = [true false];
matchTable.ColumnName = {'Selected', 'Parameter'};

applyMatch = uicontrol(fig, 'Style', 'pushbutton','String','Apply','Position',[10 10 80 30]);
applyMatch.Callback = {@applyMatching_callback,obj, sourceModel,matchTable};

hStack.sourceSMLMFit = sourceSMLMFit;
hStack.sourceModel = sourceModel;
hStack.matchTable = matchTable;
hStack.applyMatch = applyMatch;
end

function SMLMFitSelect_callback(a,b,obj,sourceModel)
selectedFitGUI = a.String{a.Value};
nameAllLoadedEval = obj.locData.SE.processors.eval.guihandles.modules.Data(:,2);
indSelected = ismember(nameAllLoadedEval, selectedFitGUI);
if any(indSelected)
obj.setPar('matchPar_fitGUIselected',indSelected);
numOfModel_selectedFitGUI = obj.locData.SE.processors.eval.processors{indSelected}.fitter.numOfModel;
% export the options
    sourceModel.String = cellfun(@num2str, num2cell(1:numOfModel_selectedFitGUI), 'UniformOutput',false);
else
    sourceModel.String = {'--'};
end
end

function parMatcher_callback(a,b,obj)
htable = obj.guihandles.anchorConvert;
data = htable.Data;
selectedRow = obj.getPar(['selectedRowConvert_' obj.name]);
fitter = obj.fitter;
try
    data(selectedRow,:) = [];
catch
end
htable.Data = data;
obj.guihandles.anchorConvert=htable;
end

function applyMatching_callback(a,b,obj, sourceModel,matchTable)
% Match the parameters with the same names.
tabConvert = obj.guihandles.anchorConvert.Data;
indSelected = obj.getPar('matchPar_fitGUIselected');
sourceSMLMModelFit = obj.locData.SE.processors.eval.processors{indSelected}.fitter;
lModel = sourceSMLMModelFit.allParsArg.model == sourceModel.Value;
parType = sourceSMLMModelFit.allParsArg.type(lModel);
parName = sourceSMLMModelFit.allParsArg.name(lModel);
parID_source = join([parType,parName],'.');
lTargetSelected = [matchTable.Data{:,1}];
targetSelected = matchTable.Data(find(lTargetSelected),2);
targetSelected_shourt = regexprep(targetSelected,'m\d\.','');
[lFound,idxSource] = ismember(targetSelected_shourt,parID_source);
ID_target = targetSelected(lFound);
ID_source = parID_source(idxSource(lFound));
prefix = ['pars.m',num2str(sourceModel.Value)];
ID_source = join([cellstr(repmat(prefix,size(ID_source))), ID_source], '.');
fitGUI_source = obj.locData.SE.processors.eval.guihandles.modules.Data{find(indSelected),2};
sizeOfNewRules = size(ID_source);
if isempty(tabConvert)
    tabConvert = [cellstr(repmat(fitGUI_source,sizeOfNewRules)), ID_source, ID_target, cell(sizeOfNewRules)];
else
    tabConvert(end+1:end+sizeOfNewRules(1),1) = cellstr(repmat(fitGUI_source,sizeOfNewRules));
    tabConvert(end-sizeOfNewRules(1)+1:end,2) = ID_source;
    tabConvert(end-sizeOfNewRules(1)+1:end,3) = ID_target;
    tabConvert(end-sizeOfNewRules(1)+1:end,4) = cell(sizeOfNewRules);
end
obj.guihandles.anchorConvert.Data = tabConvert;
end