classdef SMLMModelFit_saveResult<interfaces.DialogProcessor&interfaces.SEProcessor
% Export the results of the LocMoFitGUI
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
            obj.setPar('LocMoFit_saveResult',obj)
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
            %% This function loads data including LocMoFitter and sites.
            % basic info.
            se = obj.locData.SE;
            sites = se.sites;
            if isempty(sites)  %Jonas: this plugin stops loading SMAP in Linux
                return
            end
            roiSize = se.P.par.se_siteroi.content;
            
            listOfModules = se.processors.eval.guihandles.modules.Data(:,2);
            lModuleToSave = startsWith(listOfModules, 'LocMoFitGUI');
            nameModuleToSave = listOfModules(lModuleToSave);
            
            firstSite = sites(1);
            if isfield(firstSite.evaluation, 'generalStatistics')
                nameModuleToSave = [{'generalStatistics'}; nameModuleToSave];
                pre = 1;    % number of the modules before the LocMoFit.
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
            LocMoFitGUI_obj = se.processors.eval.processors(lModuleToSave);
            for l = 1:length(nameModuleToSave)-pre
                output.LocMoFit.(nameModuleToSave{pre+l}) = LocMoFitGUI_obj{l}.fitter;
            end
            
            % Export
            objName = class(obj);
            objName = split(objName, '.');
            
            warning('Here for now use CME3D_manager to replace LocMoFit_manager.')
            obj.fit_manager = CME3D_manager(output.sites, output.LocMoFit);
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
            allLocMoFitGUI = startsWith(namesAllEval,'LocMoFitGUI');
            names = namesAllEval(allLocMoFitGUI);
        end
        
        function registerSites_callBack(obj,a,b)
            answer = questdlg('registerSites will change the use annotation and re-sort your sites. Do you want to continue?', ...
                'Register sites', ...
                'Yes','No','No');
            % Handle response
            switch answer
                case 'Yes'
                case 'No'
                    disp('Stopped by the user.')
                    return
            end
            
            p = obj.getGuiParameters;
            sites = obj.SE.sites;
            if any([p.onlyPositive p.withoutClouds])
                fn = fieldnames(obj.locData.SE.processors.eval.children);
                lFitterGUI = strcmp('LocMoFitGUI_2',fn);
                fitter = obj.locData.SE.processors.eval.processors{lFitterGUI}.fitter;
                [~,idxCurvature] = fitter.getVariable('m1.curvature');
%                 idx = 
                for k = obj.SE.numberOfSites:-1:1
                    path = ['sites(k).evaluation.LocMoFitGUI_2' idxCurvature];
                    curvature(k) = eval(path);
                end
            end
                        
            % enable the first two sorts
            g = obj.getPar('mainGui');
            sortROIs = g.children.guiSites.children.Helper.children.SortROIs;
            sortROIs.guihandles.direction1.Value = 2;
            sortROIs.guihandles.sortedit1.String = 'annotation.use';
            sortROIs.guihandles.sortprop1.Value = 3;

            sortROIs.guihandles.direction2.Value = 1;
            sortROIs.guihandles.sortprop2.Value = 4;
            sortROIs.guihandles.sortedit2.String = 'evaluation.LocMoFitGUI_2.allParsArg.value(13)';

            % disable all other sorts
            sortROIs.guihandles.sortedit3.String = '';
            sortROIs.guihandles.sortprop3.Value = 1;
            sortROIs.guihandles.sortedit4.String = '';
            sortROIs.guihandles.sortprop4.Value = 1;

            sortROIs.run(sortROIs.getAllParameters);
            
            obj.loadData;
            obj.fit_manager.masterAvg;
            obj.locData.regroup;
            obj.locData.filter;
        end
        
        function module = dynamicRec_callBack(obj,a,b)
            % hack an evaluate plug-in in order to use the obj.getLocs(...)
            module=plugin('ROIManager','Analyze','SMLMModelFit_dynamicRec_mCME');
            p.Vrim=100;

            module.handle=figure('MenuBar','none','Toolbar','none','Name','SMLMModelFit_dynamicRec_mCME');
            module.attachPar(obj.P);
            module.attachLocData(obj.locData);

            p.Xrim=10;
            module.setGuiAppearence(p)
            module.makeGui;

            module.linkedManager = obj.fit_manager;
        %     fdcal=figure(233);
        %     dcal=plugin('ROIManager','Analyze','LocMoFit_dynamicRec_mCME',fdcal,obj.P);
        %     dcal.attachLocData(obj.SE.locData);
        %     dcal.makeGui;

        %     obj.loadData;
        %     obj.fit_manager.dynamicRec;
        %     obj.locData.regroup;
        %     obj.locData.filter;
        end
        
        function mkMovie_callBack(obj,a,b)
            obj.loadData;
            p = obj.getAllParameters;
            [file,path] = uiputfile('*.tif', 'Save as', '');

            if file~=0
                boundCurvature = [-inf inf]; % take all sites
                if p.withoutClouds
                    boundCurvature(2) = p.cloudsThreshold;
                end
                % site filtering on sites
                lFiter = obj.fit_manager.filter('LocMoFitGUI_2.m1.curvature', boundCurvature);
                lFiter = lFiter&obj.fit_manager.useSites;
                obj.fit_manager.usedSites = lFiter;
                
                % get the number of neg. curvature sites
                lFiter = obj.fit_manager.filter('LocMoFitGUI_2.m1.curvature', [-inf 0]);
                lFiter = lFiter&obj.fit_manager.useSites;
                numOfNegCur = sum(lFiter); % this for calculating the peusdotime
                
                obj.fit_manager.mkMovie('saveTo', [path file],'numberOfSitesWithNegCur',numOfNegCur);
            else
                warning('Please specify where to save.')
            end
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

pard.plotUseOnly.object=struct('Style','checkbox','Value',1,'String','Plot use only');
pard.plotUseOnly.position=[1,3.5];
pard.plotUseOnly.Width=0.7;

% pard.tPlotUseOnly.object=struct('Style','text','String','Plot use only');
% pard.tPlotUseOnly.position=[1,3.7];
% pard.tPlotUseOnly.Width=1;



pard.registerSites.object=struct('Style','pushbutton','String','Register sites','Callback', {{@obj.registerSites_callBack}});
pard.registerSites.position=[3,3.7];
pard.registerSites.Width=1;

pard.recSites.object=struct('Style','pushbutton','String','Reconstruction','Callback', {{@obj.dynamicRec_callBack}});
pard.recSites.position=[4,3.7];
pard.recSites.Width=1;

pard.mkMovie.object=struct('Style','pushbutton','String','Make movie','Callback', {{@obj.mkMovie_callBack}});
pard.mkMovie.position=[6,3.7];
pard.mkMovie.Width=1;

pard.t_winSize.object=struct('Style','text','String','Window size');
pard.t_winSize.position=[7,3.7];
pard.t_winSize.Width=0.7;

pard.winSize.object=struct('Style','edit','String','30');
pard.winSize.position=[7,4.4];
pard.winSize.Width=0.5;

pard.t_stepSize.object=struct('Style','text','String','Step size');
pard.t_stepSize.position=[8,3.7];
pard.t_stepSize.Width=0.7;

pard.onlyPositive.object=struct('Style','checkbox','Value',0,'String','Only curvature>0');
pard.onlyPositive.position=[9,3.7];
pard.onlyPositive.Width=1;

pard.withoutClouds.object=struct('Style','checkbox','Value',0,'String','No clouds');
pard.withoutClouds.position=[10,3.7];
pard.withoutClouds.Width=1;

pard.cloudsThreshold.object=struct('Style','edit','String',0.016);
pard.cloudsThreshold.position=[10,4.3];
pard.cloudsThreshold.Width=0.5;

pard.stepSize.object=struct('Style','edit','String','30');
pard.stepSize.position=[8,4.4];
pard.stepSize.Width=0.5;

pard.parsTable.object=struct('Style','text','String','table pos');
pard.parsTable.position=[12,1];
pard.parsTable.Width=2.5;
pard.parsTable.Height=10;

% pard.syncParameters={{'roimanager_processors','parsTable',{'Data'}}};

pard.plugininfo.description='Manage and export results of LocMoFitGUI.';
pard.plugininfo.type='ROI_Analyze';
end

function load_callBack(a,b,obj)
    obj.loadData;
end

function extLoad_callBack(a,b,obj)
    obj.loadData;
end





function recSettings_callBack(a,b,obj)
    fig = figure;
    set(fig,'Tag', 'recSettings', 'Name', 'Settings - Dyanmic reconstruction')
    hOld = uicontrol(fig,'Position',[50 100 100 200]); 
    fitter = obj.fit_manager.data.fitter.LocMoFitGUI_2;
    fitter.createConvertTable(hOld, 'hTable');
    
    uicontrol(fig,'Position',[50 70 25 25], 'Style', 'pushbutton','String','+','Callback',{@fitter.addRow, 'hTable'}); 
    uicontrol(fig,'Position',[75 70 25 25], 'Style', 'pushbutton','String','-','Callback',{@fitter.rmRow, 'hTable'});
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
