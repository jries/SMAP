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
                if isempty(val)
                    continue
                end
                val(isinf(val)|isnan(val))=[];
                
                switch p.setrange.Value
                    case 1
                        q=[min(val) max(val)];
                    case 2
                        q=quantile(val,[p.quantile(1), 1-p.quantile(end)]);
                    case 3     
                        q=p.quantile; 
                        if length(q)==1
                            q=[0 q];
                        end
                end
                
                if p.setbinwidth
                    BinWidth=p.binwidth;
                else
                    BinWidth=(q(2)-q(1))/200;
                end
                n=(q(1)-BinWidth:BinWidth:q(2)+BinWidth)';
                if length(n)>1e5
                    warndlg('histogram too large, set range or dont set binwidth')
                    return
                end
                h=histogram(axis1,val,n,'DisplayStyle',p.plotformat.selection);
                hold(axis1,'on')
                legends{2*k-1}=['Layer ' num2str(layers(k))];
                inval=val>=q(1) & val<=q(end);
                
                %modal value
                [~, indmax]=max(h.Values);
                fitr=2; range=(max(1,indmax-fitr):min(length(h.Values),indmax+fitr))';
                nrange=n(range)+BinWidth/2;hrange=h.Values(range)';
                fitp=fit(nrange,hrange,'poly2');
                maxval=-fitp.p2/fitp.p1/2;    
                plot(axis1,nrange,fitp(nrange))
                legends{2*k}=['fit ' num2str(maxval,2)];
                out(:,k)=[median(val(inval)),mean(val(inval)),std(val(inval)),maxval];
            end
            
                

           legend(axis1,legends);
           xlim([n(2) n(end)]);
           
           h=uitable(ax2.Parent, 'Data',out,'ColumnName',legends,'RowName',{'median','mean','std','modal'},'Units','normalized','Position',[0 0 1 1]);
           tab=ax2.Parent;
           delete(ax2)
           tab.Parent.Parent.Renderer='painters';
            
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
            pard.binwidth.position=[1,2.9];
            pard.binwidth.Width=.3;
            
            pard.setrange.object=struct('String',{{'automatic range','range (quantile)','range (absolute)'}},'Style','popupmenu','Value',2);
            pard.setrange.position=[1,3.2];
            pard.setrange.Width=1.3;
            pard.quantile.object=struct('String','0.001','Style','edit');
            pard.quantile.position=[1,4.5];
            pard.quantile.Width=0.5;
            
            pard.plotformat.object=struct('String',{{'bar','stairs'}},'Style','popupmenu');
            pard.plotformat.position=[2,1];
            pard.plotformat.Width=1.5;
            
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