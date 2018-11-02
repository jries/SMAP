classdef BatchAnalyze<interfaces.DialogProcessor&interfaces.SEProcessor
    properties
    end
    methods
        function obj=BatchAnalyze(varargin)        
            obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={};
        end
        function initGui(obj)
%             g=obj.getPar('mainGui');
            plugins=obj.getPar('menu_plugins');
            fn=fieldnames(plugins.ROIManager.Analyze);
            obj.guihandles.analyzer.String=fn;
        end
        
        function out=run(obj,p)  
            out=[];
            g=obj.getPar('mainGui');
            
            [f,pfad]=uigetfile('*.mat','MultiSelect','on');
            if ~iscell(f)
                f={f};
            end
            rp=obj.getPar('ROI_restorparamters');
            obj.setPar('ROI_restorparamters',false)
            for files=1:length(f)
                gf=g.children.guiFile;
                
                gf.loadbutton_callback(0,0,0,pfad,f{files});
                if p.redrawall

                    g.locData.SE.processors.preview.redrawall(true); %%% remove (true) to also redraw files and cells. feb 3 2018 /mm

                end
                analzyers=fieldnames(g.children.guiSites.children.Analyze.children);
                selectedan=p.analyzer.selection;
                if any(strcmp(analzyers,selectedan))
                    analyzering=g.children.guiSites.children.Analyze.children.(selectedan);
                    analyzering.processgo;
                else
                    disp('selected analyzer should be added to ROImanager/Analyzers')
                end

                [~,filen]=fileparts([pfad  f{files}]);
                if length(filen)>20
                    folder=filen(1:20);
                else
                    folder=filen;
                end
                
                if length(f{files})>60
                    fileprefix=f{files}(1:60);
                else
                    [~,fileprefix]=fileparts(f{files});
                end
                mkdir([pfad  folder]);
                g.locData.savelocs([pfad  folder filesep 'rois_' fileprefix '_sml.mat']);

               
                analyzering.saveall([pfad  folder filesep],fileprefix);
            end
            obj.setPar('ROI_restorparamters',rp);
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end


function pard=guidef(obj)


pard.t1.object=struct('String','select analyzer','Style','text');
pard.t1.position=[1,1];
pard.analyzer.object=struct('String','empty','Style','popupmenu');
pard.analyzer.position=[1,2];
pard.analyzer.Width=3;

pard.redrawall.object=struct('String','redraw all','Style','checkbox','Value',1);
pard.redrawall.position=[2,1];

pard.plugininfo.type='ROI_Analyze';
end