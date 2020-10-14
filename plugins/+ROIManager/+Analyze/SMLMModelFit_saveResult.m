classdef SMLMModelFit_saveResult<interfaces.DialogProcessor&interfaces.SEProcessor
% Export the results of the SMLMModelFitGUI
    properties
        fit_manager
        variableTable_handle
    end
    properties (Dependent)
        variableIDs
    end
    methods
        function obj=SMLMModelFit_saveResult(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'se_viewer'};
            obj.showresults=true;
        end
        
        function makeGui(obj,varargin)
            makeGui@interfaces.DialogProcessor(obj); %make the main GUI from the guidef definitions
                %Settings               
            
            obj.loadData
        end
        
        function out=run(obj,p)
            % ask user to specify the file name
            % Todo: the UpdateFcn only accept one set of x y data now. This
            % overwrites the info about earlier plots. This has to be
            % fixed.
            
            hList = obj.fit_manager.plot;
            
            fn = fieldnames(hList);
            for k = 1:length(fn)
                oldParent = hList.(fn{k}).Parent;
                dc = datacursormode(oldParent);
                
                tab = obj.initaxis(fn{k});
                newParent = tab.Parent;
                dcNew = datacursormode(newParent.Parent.Parent);
                
                dcNew.UpdateFcn = dc.UpdateFcn;
                delete(tab)
                hList.(fn{k}).Parent = newParent;
                delete(oldParent)
            end
                        
            out = [];
        end
        
        function loadData(obj,p)
            %% This function loads data including SMLMModelFitter and sites.
            % basic info.
            se = obj.locData.SE;
            sites = se.sites;
            roiSize = se.P.par.se_siteroi.content;
            
            listOfModules = se.processors.eval.guihandles.modules.Data(:,2);
            lModuleToSave = startsWith(listOfModules, 'SMLMModelFitGUI');
            nameModuleToSave = listOfModules(lModuleToSave);
            
            firstSite = sites(1);
            if isfield(firstSite.evaluation, 'generalStatistics')
                nameModuleToSave = [{'generalStatistics'}; nameModuleToSave];
                pre = 1;    % number of the modules before the SMLMModelFit.
            else
                pre = 0;
            end
                
            export_sites = [];
            
            % Site info
            for k = se.numberOfSites:-1:1
                oneSite = sites(k);

                export_sites(k).annotation.use = oneSite.annotation.use;
                export_sites(k).pos = oneSite.pos;
                export_sites(k).ID = oneSite.ID;
                export_sites(k).info = oneSite.info;
                for l = 1:length(nameModuleToSave)
                    export_sites(k).evaluation.(nameModuleToSave{l}) = oneSite.evaluation.(nameModuleToSave{l});
                end
            end
            sites = export_sites;
            output.sites = sites;
            
            
            
            % Fit info.
            SMLMModelFitGUI_obj = se.processors.eval.processors(lModuleToSave);
            for l = 1:length(nameModuleToSave)-pre
                output.SMLMModelFit.(nameModuleToSave{pre+l}) = SMLMModelFitGUI_obj{l}.fitter;
            end
            
            % Export
            objName = class(obj);
            objName = split(objName, '.');
            
            warning('Here for now use CME3D_manager to replace SMLMModelFit_manager.')
            obj.fit_manager = CME3D_manager(output.sites, output.SMLMModelFit);
            obj.fit_manager.parentObj = obj;
            
            % copy the plot setting from the exsiting fitManager
            old_manager = obj.getPar([objName{end} '_fitManager']);
            
            
            hTable =obj.guihandles.parsTable;
            IDs = obj.fit_manager.variableIDs;
            
            grps = cell(size(IDs,1),1);
            orderVariable = cell(size(IDs,1),1);
            if ~isempty(old_manager)&&~isempty(old_manager.plotSettings)
                obj.fit_manager.plotSettings = old_manager.plotSettings;
                grpsOld = plotSettings2grp(obj.fit_manager.plotSettings);
                [~,ind] = ismember(grpsOld(:,1),IDs);
                grps(ind) = grpsOld(:,2);
            end
            if ~isempty(old_manager)&&~isempty(old_manager.variableTableCol)
                [~,orderNew] = ismember(old_manager.variableTableCol,IDs);
                orderVariable(orderNew) = num2cell(1:length(orderNew));
            end
            lSelect = false(size(IDs));
            obj.variableTable_handle = uitable(hTable.Parent,'Data',[IDs grps orderVariable],'Position',hTable.Position,'ColumnWidth', {230,20,20}, 'ColumnEditable', [false true true],'CellEditCallback',{@variableTableEditCallback, obj});
            obj.variableTable_handle.ColumnName = {'Variables','Plot','Export order'};
            obj.variableTable_handle.ColumnFormat = {'char','char','char'};
            
            obj.setPar([objName{end} '_fitManager'], obj.fit_manager)
        end
              
        function plotSettings = getPlotSettings(obj)
            % restructure the plot settings.
            % plotSettings is a struct with fields of plot groups, having
            % the field names defined as the plot marker.
            grp = obj.variableTable_handle.Data(:,2);
            IDs = obj.variableTable_handle.Data(:,1);
            indGrp = find(~cellfun('isempty',grp));
            
            for k = 1:length(indGrp)
                oneGrps = split(grp{indGrp(k)}, ' '); % on variable could appear in multiple plots
                for l = 1:length(oneGrps)
                    grpParts = split(oneGrps{l}, '.');
                    oneID = IDs{indGrp(k)};
                    if length(grpParts)==2
                        switch grpParts{2}
                            case 'x'
                                indAxis = 1;
                            case 'y'
                                indAxis = 2;
                        end
                        plotSettings.(grpParts{1}){indAxis} = oneID;
                    else
                        indAxis = 1;
                        plotSettings.(grpParts{1}){indAxis} = oneID;
                    end
                end
            end
        end
        
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        
        function names = getFitterNames(obj)
            namesAllEval = fieldnames(obj.locData.SE.processors.eval.children);
            allSMLMModelFitGUI = startsWith(namesAllEval,'SMLMModelFitGUI');
            names = namesAllEval(allSMLMModelFitGUI);
        end
        
        function createPlot(obj)
        end
    end
end




function pard=guidef(obj)

pard.load.object=struct('Style','pushbutton','String','Refresh','Callback', {{@load_callBack,obj}});
pard.load.position=[1,1];
pard.load.Width=1;

pard.save.object=struct('Style','pushbutton','String','Save','Callback', {{@save_callBack,obj}});
pard.save.position=[1,2];
pard.save.Width=1;

pard.extLoad.object=struct('Style','pushbutton','String','Load from external','Callback', {{@extLoad_callBack,obj}});
pard.extLoad.position=[2,1];
pard.extLoad.Width=1;

pard.plotUseOnly.object=struct('Style','checkbox','Value',1);
pard.plotUseOnly.position=[1,3.5];
pard.plotUseOnly.Width=0.2;

pard.tPlotUseOnly.object=struct('Style','text','String','Plot use only');
pard.tPlotUseOnly.position=[1,3.7];
pard.tPlotUseOnly.Width=1;

pard.registerSites.object=struct('Style','pushbutton','String','Register sites','Callback', {{@registerSites_callBack,obj}});
pard.registerSites.position=[3,3.7];
pard.registerSites.Width=1;

pard.recSites.object=struct('Style','pushbutton','String','Reconstruction','Callback', {{@dynamicRec_callBack,obj}});
pard.recSites.position=[4,3.7];
pard.recSites.Width=1;

pard.parsTable.object=struct('Style','text','String','table pos');
pard.parsTable.position=[12,1];
pard.parsTable.Width=2.5;
pard.parsTable.Height=10;

% pard.syncParameters={{'roimanager_processors','parsTable',{'Data'}}};

pard.plugininfo.description='Manage and export results of SMLMModelFitGUI.';
pard.plugininfo.type='ROI_Analyze';
end

function load_callBack(a,b,obj)
    obj.loadData;
end

function extLoad_callBack(a,b,obj)
    obj.loadData;
end

function registerSites_callBack(a,b,obj)
    obj.loadData;
    obj.fit_manager.masterAvg;
    obj.locData.regroup;
    obj.locData.filter;
end

function dynamicRec_callBack(a,b,obj)
    obj.loadData;
    obj.fit_manager.dynamicRec;
    obj.locData.regroup;
    obj.locData.filter;
end

function variableTableEditCallback(a,b,obj)
    if b.Indices(2) == 2
        plotSettings = obj.getPlotSettings;
        obj.fit_manager.plotSettings = plotSettings;
    elseif b.Indices(2) == 3
        dataTable = obj.variableTable_handle.Data;
        
        variableTableCol = {};
        ind = str2double(dataTable(:,3));
        variableTableCol(ind(ind>0)) = dataTable(ind>0,1);
        
        obj.fit_manager.variableTableCol = variableTableCol;
        
        objName = class(obj);
        objName = split(objName, '.');
        obj.setPar([objName{end} '_fitManager'], obj.fit_manager)
    end
end

function save_callBack(a,b,obj)
    [file,path] = uiputfile('*_fitResult.mat', 'Save as', '');
    if file~=0
        obj.loadData;
        % only when the path is specified
        fit_manager = obj.fit_manager;
        file = strsplit(file,'.');
        if ~endsWith(file{1}, '_fitResult')
            fileMat = strcat(path, file{1}, '_fitResult.', file{2});
        else
            fileMat = strcat(path, file{1}, '.', file{2});
        end
        
        % save the fit_manger first
%         save(fileMat, 'fit_manager')
        
        % get which parameter(s) to save from the table
        dataTable = obj.variableTable_handle.Data;
        
        variableTableCol = {};
        if isa([dataTable{:,3}],'double')
            lEmpty = cellfun(@isempty,dataTable(:,3));
            ind = dataTable(:,3);
            ind(lEmpty) = {0};
            ind = [ind{:}];
        else
            ind = str2double(dataTable(:,3));
        end
        variableTableCol(ind(ind>0)) = dataTable(ind>0,1);
        
        obj.fit_manager.variableTableCol = variableTableCol;
        variableTable = obj.fit_manager.variableTable;
        
        objName = class(obj);
        objName = split(objName, '.');
        obj.setPar([objName{end} '_fitManager'], obj.fit_manager)
        writematrix(variableTable,strcat(path, [file{1} '.txt']),'Delimiter','tab')
    else
        warning('Please specify where to save.')
    end
end

function grp = plotSettings2grp(plotSettings)
    fn = fieldnames(plotSettings);
    grp = {};
    for k = 1:length(fn)
        IDs = plotSettings.(fn{k});
        oneGrp = {};
        for l = 1:length(IDs)
            if length(IDs)==2
                if l == 2
                    oneGrp(l,:) = {IDs{l}, [fn{k} '.y']};
                else
                    oneGrp(l,:) = {IDs{l}, [fn{k} '.x']};
                end
            else
                oneGrp(l,:) = {IDs{l}, fn{k}};
            end
        end
        grp = [grp; oneGrp];
    end
    [grpUni, ~, ind] = unique(grp(:,1));
    for k = 1:length(grpUni)
        l = ind == k;
        toBeCombined = grp(l,2);
        if length(toBeCombined)>1
            toBeCombined(:,2) = {repelem(" ",length(toBeCombined)-1) []};
            toBeCombined = toBeCombined';
            grpUni(k,2) = cellstr(strcat(toBeCombined{:}));
        else
            grpUni(k,2) = toBeCombined;
        end
    end
    grp = grpUni;
end
