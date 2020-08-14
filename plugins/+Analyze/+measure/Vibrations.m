classdef Vibrations<interfaces.DialogProcessor
    % Vibrations calculates the fourier transform of x, y and z to
    % visualize vibrations
    properties
    end
    methods
        function obj=Vibrations(varargin)   %replace by filename        
            obj@interfaces.DialogProcessor(varargin{:}) 
            obj.showresults=true; %set true, if results are shown by default
            obj.history=false; %if set true, every time the plugin is called, its parameters are added to the history. Makes sense only if plugin changes the data
        end       
        function initGui(obj)
%             obj.guihandles.locfield=obj.getPar('locFields');
        end
        function out=run(obj,p)
            layers=find(obj.getPar('sr_layerson'));
            plotfields={'xnm','ynm','znm'};
            plotfieldserr={'xnmerr','ynmerr','locprecznm'};
            for k=1:length(plotfields)
                axx(k)=obj.initaxis(plotfields{k}(1));
                axfx(k)=obj.initaxis(['fft(' plotfields{k}(1) ')']);
            end
            
            legends={};
            for k=1:length(layers)
                legends{k}=['layer' num2str(layers(k))];
                locs=obj.locData.getloc({'xnm','ynm','znm','frame','filenumber','xnmerr','ynmerr','locprecznm','time'},'position','roi','layer',layers(k),'grouping','ungrouped');
                if ~isempty(locs.time)
                    locsnew=locs;
                    told=locs.time;
                    dt=quantile(diff(told),0.1);
                    t=(min(told):dt:max(told))';
                    Fs=1000/dt;
                    for f=1:length(plotfields)
                        x=locs.(plotfields{f});
                        if ~isempty(x)
                        xn=interp1(told,x,t);
                        else
                            xn=[];
                        end
                        locsnew.(plotfields{f})=xn;
                    end
                    locs=locsnew;
                else
                    if p.overwriteframetime
                        dt=p.frametime;
                    else
                    dt=obj.locData.files.file(locs.filenumber(1)).info.timediff;
                    end
                    dt
                    Fs=1000/dt;
                    t=locs.frame*dt;
                end
                for f=1:length(plotfields)
                    x=locs.(plotfields{f});
                    if isempty(x)
                        continue
                    end
                    plot(axx(f),t,x-mean(x))
                    xlabel(axx(f),'time(ms)')
                    ylabel(axx(f),['d' plotfields{f} '(nm)']);
                    title(axx(f),['SD = ' num2str(std(x),3) ', mean(locprec) = ' num2str(mean(locs.(plotfieldserr{f})),3) ' nm']);
                    L=length(x);
                    xf=abs(fft(x)/L);
                    xfp=2*xf(1:L/2+1);
                    xfp(1)=0;
                    freq=Fs*(0:(L/2))/L;
                    semilogy(axfx(f),freq,xfp);
                    xlabel(axfx(f),'frequency (Hz)')
                    ylabel(axfx(f),'Amplitude (nm), 2*fft (x)');
                    
                    hold(axx(f),'on')
                    hold(axfx(f),'on')
                    [m,ind]=max(xfp);
                    mfreq=freq(ind);
                    title(axfx(f),['maximum at ' num2str(mfreq) ' Hz']);
                    plot(axfx(f),mfreq,m,'x')
                end

            end
  
            for f=1:length(plotfields)
                legend(axx(f),legends)
                legend(axfx(f),legends)
            end
    
            
            out=[]; %no output
%             out.clipboard={'results1',3,'text1'}; % out.clipboard is copied to clipboard, separated by tabs.
%             out.error='this error occured because you did something wrong'; %if an error occoured, you can output it in the status bar with this command.
        end
        function pard=guidef(obj)
            pard.text.object=struct('String','Analyze vibration measurements. Only 1 bead in ROI.','Style','text');
            pard.text.position=[1,1];
            pard.text.Width=4;
            pard.text2.object=struct('String','Manually click Refresh in MM property browser after changing exposure time.','Style','text');
            pard.text2.position=[2,1];
            pard.text2.Width=4;
            pard.overwriteframetime.object=struct('String','Overwrite frame time, set to (ms):','Style','checkbox');
            pard.overwriteframetime.position=[3,1];
            pard.overwriteframetime.Width=2;
            pard.frametime.object=struct('String','20','Style','edit');
            pard.frametime.position=[3,3];
            pard.frametime.Width=0.5;
%             pard.setbinwidth.object=struct('String','set binwidth','Style','checkbox','Value',0);
%             pard.setbinwidth.position=[1,2];
%             pard.binwidth.object=struct('String','','Style','edit','Value',0);
%             pard.binwidth.position=[1,3];
%             pard.binwidth.Width=.5;
%             
%             pard.setrange.object=struct('String','range (quantile)','Style','checkbox');
%             pard.setrange.position=[1,3.5];
%             pard.quantile.object=struct('String','0.001','Style','edit');
%             pard.quantile.position=[1,4.5];
%             pard.quantile.Width=0.5;
            pard.plugininfo.name='Vibrations';
            pard.plugininfo.description='Calculates the fourier transform of x, y and z to visualize vibrations';
            pard.plugininfo.type='ProcessorPlugin'; %type of plugin. Currently: ProcessorPlugin, WorkflowModule, WorkflowFitter, Renderer, LoaderPlugin, SaverPlugin, ROI_Analyze, ROI_Evaluate,WorkflowIntensity
  
%             pard.syncParameters={{'locFields','locfield',{'String'}}};
            

        end
    end
end


function callbackfunction(uiobject,data,obj,extradata)
disp('callback')
end

function aftersync_callback(obj)
disp('aftersync')
end