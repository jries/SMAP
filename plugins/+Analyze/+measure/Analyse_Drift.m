classdef Analyse_Drift<interfaces.DialogProcessor
    % Analyse_Drift analyses and displays drift from time-lapse images of beads
    properties
    end
    methods
        function obj=Analyse_Drift(varargin)   %replace by filename        
            obj@interfaces.DialogProcessor(varargin{:}) 
            obj.showresults=true; %set true, if results are shown by default
            obj.history=false; %if set true, every time the plugin is called, its parameters are added to the history. Makes sense only if plugin changes the data
        end       
        function initGui(obj)
%             obj.guihandles.locfield=obj.getPar('locFields');
        end
        function out=run(obj,p)

            % identify beads
            % combine split beads into single beads
            % plot all beads vs time
            % calculate meaningful average: even if beads don't persist for
            % entire time.
            % Plot average, maybe add smoothed average
       
            layers=find(obj.getPar('sr_layerson'));

            locs=obj.locData.getloc({'znm'},'position','roi','layer',1,'grouping','ungrouped');
            if isempty(locs.znm)
                plotfields={'xnm','ynm'};
            else
                plotfields={'xnm','ynm','znm'};
            end
        

                locs=obj.locData.getloc({'phot','xnm','ynm','znm','locprecxnm','locprecznm','frame','filenumber','time','numberInGroup','groupindex'},'position','roi','layer',layers,'grouping','ungrouped');
                beads=unique(locs.groupindex(locs.numberInGroup>p.minframes));
                maxframes=max(locs.frame);
                framesall=(1:maxframes)';
                numbeads=length(beads);

                if p.overwriteframetime
                    dt=p.frametime;
                else
                    dt=obj.locData.files.file(locs.filenumber(1)).info.timediff;
                end

                dxall=zeros(maxframes,numbeads,3);
                norm=zeros(maxframes,numbeads,3);
                dxav=zeros(maxframes,1,3);
                phot=zeros(maxframes,numbeads);


                for b=1:length(beads)
                    thisbead=locs.groupindex==beads(b);

                    xh=locs.xnm(thisbead);yh=locs.ynm(thisbead);
                    if isempty(locs.znm)
                        zh=0*xh;
                    else
                        zh=locs.znm(thisbead);
                    end
                    coordsh=[];
                    coordsh(:,1,:)=[xh,yh,zh];
                    
                    frameh=locs.frame(thisbead);
                    framegood=find(sum(norm(:,:,1),2)>0); %look only at x, 
                    if isempty(framegood) %first bead
                        framegood=frameh;
                    end
                    [frameref, inh]=intersect(frameh,framegood);

                   %get reference averages
                    meanref=mean(dxav(frameref,1,:),1);
                    avxh=mean(coordsh(inh,1,:),1);
                    dxh=coordsh-avxh+meanref;
                    dxall(frameh,b,:)=dxh;
                    norm(frameh,b,:)=1;
                    dxav=sum(dxall.*norm,2)./sum(norm,2);
                    phot(frameh,b)=locs.phot(thisbead);

                    
                end
                minframe=find(sum(norm(:,:,1),2)>0,1,'first');
                timeplot=(framesall-minframe)*dt;

                axall=obj.initaxis('average');

                for pf=1:length(plotfields)
                    axx(pf)=obj.initaxis(plotfields{pf}(1));
                    hold(axx(pf),'off')
                    for b=1:numbeads
                        isgood=norm(:,b,pf)>0;
                        plot(axx(pf),timeplot(isgood),dxall(isgood,b,pf))
                        hold(axx(pf),'on')
                    end
                    plot(axx(pf),timeplot,dxav(:,:,pf),'k','LineWidth',2);
                    xlabel('time (s)')
                    ylabel('drift (nm)')

                    plot(axall,timeplot,dxav(:,:,pf));
                    hold(axall,'on')

                    
                end
                legend(axall,plotfields)
                    xlabel(axall,'time (s)')
                    ylabel(axall,'drift (nm)')

                    axxp=obj.initaxis('phot');
                    hold(axxp,'off')
                    for b=1:numbeads
                        % isgood=norm(:,b,pf)>0;
                        plot(axxp,timeplot,phot(:,b))
                        hold(axxp,'on')
                    end
                    % plot(axxp,timeplot,photav(:,:,pf),'k','LineWidth',2);
                    xlabel('time (s)')
                    ylabel('drift (nm)')

                    plot(axall,timeplot,dxav(:,:,pf));
                    hold(axall,'on')
%                 if ~isempty(locs.time) %MINFLUX or similar
%                     locsplot=locs;
%                     locsnew=locs;
%                     told=locs.time;
%                     tplot=told;
%                     dt=quantile(diff(told),0.2)/2
%                     t=(min(told):dt:max(told))';
%                     Fs=1000/dt;
%                     for f=1:length(plotfields)
%                         x=locs.(plotfields{f});
%                         if ~isempty(x)
%                         xn=interp1(told,x,t);
%                         else
%                             xn=[];
%                         end
%                         locsnew.(plotfields{f})=xn;
%                     end
%                     locs=locsnew;
%                 else %

            out=[]; %no output
        end
        function pard=guidef(obj)
            pard.text.object=struct('String','Analyze drift from beads. ','Style','text');
            pard.text.position=[1,1];
            pard.text.Width=4;
%             pard.text2.object=struct('String','Manually click Refresh in MM property browser after changing exposure time.','Style','text');
%             pard.text2.position=[2,1];
%             pard.text2.Width=4;

            pard.overwriteframetime.object=struct('String','Overwrite frame time, set to (ms):','Style','checkbox');
            pard.overwriteframetime.position=[3,1];
            pard.overwriteframetime.Width=2;
            pard.frametime.object=struct('String','1000','Style','edit');
            pard.frametime.position=[3,3];
            pard.frametime.Width=0.5;


            pard.minframest.object=struct('String','Minimum length (frames) to be counted as bead','Style','text');
            pard.minframest.position=[4,1];
            pard.minframest.Width=2;
            pard.minframes.object=struct('String','100','Style','edit');
            pard.minframes.position=[4,3];
            pard.minframes.Width=0.5;

            pard.plugininfo.name='Drift analysis';
            pard.plugininfo.description='Calculates and visualizes drift from berads';
            pard.plugininfo.type='ProcessorPlugin'; %type of plugin. Currently: ProcessorPlugin, WorkflowModule, WorkflowFitter, Renderer, LoaderPlugin, SaverPlugin, ROI_Analyze, ROI_Evaluate,WorkflowIntensity
        end
    end
end
