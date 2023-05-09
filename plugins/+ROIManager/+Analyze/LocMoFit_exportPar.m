classdef LocMoFit_exportPar<interfaces.DialogProcessor&interfaces.SEProcessor
    % Export the results of the LocMoFitGUI
    properties
        fit_manager
        variableTable_handle
    end
    properties (Dependent)
        variableIDs
    end
    methods
        function obj=LocMoFit_exportPar(varargin)
            obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'se_viewer'};
            obj.showresults=true;
        end

        function makeGui(obj,varargin)
            makeGui@interfaces.DialogProcessor(obj); %make the main GUI from the guidef definitions
            %Settings

            obj.loadData
            obj.updateOptions
            obj.setPar('LocMoFit_exportPar',obj)
        end

        function out=run(obj,p)
            % ask user to specify the file name
            % Todo: the UpdateFcn only accept one set of x y data now. This
            % overwrites the info about earlier plots. This has to be
            % fixed.
            
            lastfile=obj.getPar('lastSMLFile');
            defaultName = [lastfile(1:end-8) '_' p.which_LocMoFitGUI.selection '_parameters.txt'];

            [file,path] = uiputfile('*.txt', 'Save as', defaultName);
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

                variableTable = obj.defineWhat2Load(p);

                objName = class(obj);
                objName = split(objName, '.');
                obj.setPar([objName{end} '_fitManager'], obj.fit_manager)
                writematrix(variableTable,strcat(path, [file{1} '.txt']),'Delimiter','tab')
            else
                warning('Please specify where to save.')
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

            obj.fit_manager = LocMoFit_manager(output.sites, output.LocMoFit);
            obj.fit_manager.parentObj = obj;

            obj.setPar([objName{end} '_fitManager'], obj.fit_manager)
        end

        function updateOptions(obj)
            names_LocMoFitGUI = obj.getFitterNames;
            if isempty(names_LocMoFitGUI)
                warning('At least one LocMoFitGUI instance has to be defined.')
            else
                obj.guihandles.which_LocMoFitGUI.String = names_LocMoFitGUI;
            end
        end

        function variableTable = defineWhat2Load(obj, p)
            switch num2str([p.l_mPar p.l_lPar p.l_dPar])
                case '1  1  1'
                    type = 'all';
                case '1  1  0'
                    type = 'main';
                case '1  0  0'
                    type = 'mPar';
                case '0  1  0'
                    type = 'lPar';
                case '0  0  1'
                    type = 'auxiliary';
            end
            fitterName = p.which_LocMoFitGUI.selection;
            variableIDs = obj.fit_manager.variableIDs('fitterNames',fitterName,'type', type);
            obj.fit_manager.variableTableCol = variableIDs';  
            variableTable = obj.fit_manager.variableTable;
        end

        function pard=guidef(obj)
            pard=guidef(obj);
        end

        function names = getFitterNames(obj)
            namesAllEval = fieldnames(obj.locData.SE.processors.eval.children);
            allLocMoFitGUI = startsWith(namesAllEval,'LocMoFitGUI');
            names = namesAllEval(allLocMoFitGUI);
        end

    end
end

function pard=guidef(obj)

pard.t_which_LocMoFitGUI.object=struct('Style','text','String','Which LocMoFitGUI?');
pard.t_which_LocMoFitGUI.position=[1,1];
pard.t_which_LocMoFitGUI.Width=1;

pard.which_LocMoFitGUI.object=struct('Style','popupmenu','Value',1,'String','Plot use only');
pard.which_LocMoFitGUI.position=[1,2];
pard.which_LocMoFitGUI.Width=0.7;

pard.t_whatToExport.object=struct('Style','text','String','What to export:');
pard.t_whatToExport.position=[2,1];
pard.t_whatToExport.Width=1;

pard.l_siteInfo.object=struct('Style','checkbox','String','Information of sites', 'Value', 1);
pard.l_siteInfo.position=[3,1];
pard.l_siteInfo.Width=1;

pard.l_mPar.object=struct('Style','checkbox','String','Intrinsic parameters', 'Value', 1);
pard.l_mPar.position=[3,2];
pard.l_mPar.Width=1;

pard.l_lPar.object=struct('Style','checkbox','String','Extrinsic parameters', 'Value', 0);
pard.l_lPar.position=[3,3];
pard.l_lPar.Width=1;

pard.l_dPar.object=struct('Style','checkbox','String','Derived parameters', 'Value', 0);
pard.l_dPar.position=[3,4];
pard.l_dPar.Width=1;

pard.plugininfo.description='Export parameter values LocMoFitGUI.';
pard.plugininfo.type='ROI_Analyze';
end

