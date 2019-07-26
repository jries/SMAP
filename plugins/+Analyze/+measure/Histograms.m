classdef Histograms<interfaces.DialogProcessor
    % Histograms calculates histograms and statistics for any localization attribute (field)  
    properties
    end
    methods
        function obj=Histograms(varargin)   %replace by filename        
            obj@interfaces.DialogProcessor(varargin{:}) 
            obj.showresults=true; %set true, if results are shown by default
            obj.history=false; %if set true, every time the plugin is called, its parameters are added to the history. Makes sense only if plugin changes the data
        end       
        function initGui(obj)
%             obj.guihandles.locfield=obj.getPar('locFields');
        end
        function out=run(obj,p)
            layers=find(obj.getPar('sr_layerson'));
            ax2=obj.initaxis('Statistics');
            axis1=obj.initaxis('histograms');
            
            for k=1:length(layers)
                locs=obj.locData.getloc({p.locfield.selection},'position','roi','layer',layers(k));
                val=locs.(p.locfield.selection);
                h=histogram(axis1,val);
                if p.setbinwidth
                    h.BinWidth=p.binwidth;
                end
                if p.setrange
                    q=quantile(val,[p.quantile, 1-p.quantile]);
                    h.BinLimits=q;
                end
                hold(axis1,'on')
                legends{k}=['Layer ' num2str(k)];
                out(:,k)=[median(val),mean(val),std(val)];
            end
           legend(axis1,legends);
           
           h=uitable(ax2.Parent, 'Data',out,'ColumnName',legends,'RowName',{'median','mean','std'},'Units','normalized','Position',[0 0 1 1]);
           delete(ax2)
            
    
            
            out=[]; %no output
            out.clipboard={'results1',3,'text1'}; % out.clipboard is copied to clipboard, separated by tabs.
%             out.error='this error occured because you did something wrong'; %if an error occoured, you can output it in the status bar with this command.
        end
        function pard=guidef(obj)
            pard.locfield.object=struct('String',{{' '}},'Style','popupmenu','Value',1);
            pard.locfield.position=[1,1];
            
            pard.setbinwidth.object=struct('String','set binwidth','Style','checkbox','Value',0);
            pard.setbinwidth.position=[1,2];
            pard.binwidth.object=struct('String','','Style','edit','Value',0);
            pard.binwidth.position=[1,3];
            pard.binwidth.Width=.5;
            
            pard.setrange.object=struct('String','range (quantile)','Style','checkbox');
            pard.setrange.position=[1,3.5];
            pard.quantile.object=struct('String','0.001','Style','edit');
            pard.quantile.position=[1,4.5];
            pard.quantile.Width=0.5;
            pard.plugininfo.name='Histograms';
            pard.plugininfo.description='Calculates histograms and statistics for any localization attribute (field)';
            pard.plugininfo.type='Histograms'; %type of plugin. Currently: ProcessorPlugin, WorkflowModule, WorkflowFitter, Renderer, LoaderPlugin, SaverPlugin, ROI_Analyze, ROI_Evaluate,WorkflowIntensity
  
            pard.syncParameters={{'locFields','locfield',{'String'}}};
            

        end
    end
end


function callbackfunction(uiobject,data,obj,extradata)
disp('callback')
end

function aftersync_callback(obj)
disp('aftersync')
end