classdef SEMainGui< gui.GuiPluginWindow
    properties
        SEpreview
%         maindir
%         plugindir
    end
    methods
        function obj=SEMainGui(varargin)
            obj@gui.GuiPluginWindow(varargin{:})
        end
        function makeGui(obj)
            global se
            h.sitetabs=uitabgroup(obj.handle);
            h.tab_settings=uitab(h.sitetabs,'Title','Settings');
            h.tab_annotate=uitab(h.sitetabs,'Title','Annotation');
%             h.tab_sort=uitab(h.sitetabs,'Title','Sort');
            h.tab_evaluate=uitab(h.sitetabs,'Title','Evaluate');
%             h.tab_analyze=uitab(h.sitetabs,'Title','Analyze');
            obj.guihandles=h;
            
            se=interfaces.SiteExplorer;
            se.attachPar(obj.P);
            se.attachLocData(obj.locData);
%             obj.attachSE(se);
            obj.locData.SE=se;

            %settings panel
            guipar.FieldHeight=22;
            h.settingspanel=uipanel(h.tab_settings,'Unit','pixels','Position',obj.guiPar.tabsize2); 
            settings=gui.SEGUISettings(h.settingspanel,obj.P);
            settings.attachSE(se);
            settings.makeGui;
            se.processors.settings=settings;
            obj.children.settings=settings;
            
            %annotation panel
            h.annotationpanel=uipanel(h.tab_annotate,'Unit','pixels','Position',obj.guiPar.tabsize2); 
            annotation=gui.SEAnnotation(h.annotationpanel,obj.P);
            annotation.attachSE(se);
            annotation.attachLocData(obj.locData);
            annotation.setGuiAppearence(guipar);
            annotation.makeGui;
            se.processors.annotation=annotation;
            obj.children.annotation=annotation;
            
            %sort panel
%             guiparsort.FieldHeight=24;
%              h.sortpanel=uipanel(h.tab_sort,'Unit','pixels','Position',obj.guiPar.tabsize2); 
%             sort=gui.SESort(h.sortpanel,obj.P);
%             sort.attachSE(se);
%             sort.attachLocData(obj.locData);
%             sort.setGuiAppearence(guiparsort);
%             sort.makeGui;
%              se.processors.sort=sort;
%             obj.children.sort=sort;
            
            %evaluation panel
            h.evalpanel=uipanel(h.tab_evaluate,'Unit','pixels','Position',obj.guiPar.tabsize2); 
            eval=gui.SEEvaluationGui(h.evalpanel,obj.P);
            eval.attachSE(se);
            eval.attachLocData(obj.locData);
            eval.setGuiAppearence(guipar);
            eval.makeGui;
             se.processors.eval=eval;
            
            obj.children.eval=eval;
            
            %Explorer GUI
            obj.make_siteexplorer;
            
            se.processors.SEMainGui=obj;
            obj.tabgroup=h.sitetabs;
             makeGui@gui.GuiPluginWindow(obj);

        end
        
        function make_siteexplorer(obj)
            if isempty(obj.SEpreview)||~isvalid(obj.SEpreview.handle)
                pfig=figure;
                pfig.Visible='off';
                delete(pfig.Children)
                SEpreview=gui.SEExploreGui(pfig,obj.P);
                obj.locData.SE.processors.preview=SEpreview;
                SEpreview.attachLocData(obj.locData);
                SEpreview.attachSE(obj.locData.SE);
                SEpreview.makeGui;
                obj.SEpreview=SEpreview;
                
                obj.setPar('se_viewer',SEpreview);
            end

        end
    end
end