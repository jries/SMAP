classdef StepsMINFLUX<interfaces.SEEvaluationProcessor
%     
    properties
        peakline
        roihandle
        axstep
    end
    methods
        function obj=StepsMINFLUX(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj, p)
           if isfield(obj.site.evaluation,obj.name) 
                out=obj.site.evaluation.(obj.name);
           else
                out.range=[];
           end
           locs=obj.getLocs({'xnm','ynm','groupindex','tid'},'layer',find(obj.getPar('sr_layerson')),'size',obj.getPar('se_siteroi')/2);
           id=mode(locs.tid);
           index=obj.locData.loc.tid==id;
           x=obj.locData.loc.xnm(index);
           y=obj.locData.loc.ynm(index);
           time=obj.locData.loc.time(index);

           %find direction         
           c = cov(x-mean(x), y-mean(y));
           [a, ev] = eig(c);
           [ev,ind] = sort(diag(ev));
           [xa, ya] = deal(a(1,ind(end)), a(2,ind(end)));
           angle = cart2pol(xa, ya);
           [xr,yr]=rotcoord(x-mean(x),y-mean(y),angle);
           indx=xr<mean(xr);
           if mean(time(indx))>mean(time(~indx))
               xr=-xr;
           end

           axxy=obj.setoutput('xy');
           plot(axxy, xr, yr)
           axis(axxy,'equal')
           xlabel(axxy,'x (nm)')
           ylabel(axxy,'y (nm)')

           ax2=obj.setoutput('steps');
           timeplot=time-min(time);
           hold(ax2,'off')
           plot(ax2,timeplot,xr)
            hold(ax2,'on')
           xlabel(ax2,'time (ms)')
           ylabel(ax2,'position (nm)')

           %stepfinder
            if isfield(out,'range') && ~isempty(out.range)
                plot(ax2,out.range(1)*[1 1],[min(xr) max(xr)],'k--')
                plot(ax2,out.range(2)*[1 1],[min(xr) max(xr)],'k--')
                indtime=timeplot>=out.range(1) & timeplot<=out.range(2);
            else
                indtime=true(size(xr));
            end
            p.fitmean=contains(p.fitmode.selection,'mean');
            steps=AutoStepfinderRies(xr(indtime),p);

            inds=[1 ;steps.properties.IndexStep ];
            mv=[steps.properties.LevelBefore(1) ;steps.properties.LevelAfter];
            tv=timeplot(indtime);
            tv=tv(inds);
           
            stairs(ax2,tv,mv,'r')


            dmv=diff(mv);
            for k=1:length(dmv)
                text(ax2,tv(k+1),mean(mv(k:k+1)),num2str(dmv(k),'%2.0f'),'FontSize',11,'Color','magenta')
            end
            obj.axstep=ax2;

            ax3=obj.setoutput('stephist');
            stepsize=steps.properties.StepSize;
            histogram(ax3,stepsize,min(stepsize):1:max(stepsize));
            out.steps=steps;

            %step time
            steptime=diff(tv);
            out.steptime=steptime;
            axstept=obj.setoutput('steptime');
            dt=round(mean(steptime)/3);
            histogram(axstept,steptime,0:dt:max(steptime))
            title(axstept,"mean step time = "+ num2str(mean(steptime),'%2.1f') + " ms");
            % from cluster MINFLUX
            filelist=obj.getPar('filelist_short');
            out.filename=filelist.String{mode(locs.filenumber)};
            dt=diff(time);
            dtmin=min(dt);
            dtmedian=median(dt);
            dtmean=mean(dt);

            x=obj.locData.loc.xnm(index);

            efo=median(obj.locData.loc.efo(index));
            cfr=median(obj.locData.loc.cfr(index));
            eco=median(obj.locData.loc.eco(index));
            ecc=median(obj.locData.loc.ecc(index));
            efc=median(obj.locData.loc.efc(index));
            nlocs=length((index));
            ontime=max(time)-min(time);

            ltime=time-min(time);

            sigmax=std(xr);sigmay=std(yr);
            sxdetrend=std(diff(xr))/sqrt(2);sydetrend=std(diff(yr))/sqrt(2);
            [~, sxrobust]=robustMean(locs.xnm); [~, syrobust]=robustMean(locs.ynm);
            
            %graphs
            ff='%2.1f';
%             axx=obj.initaxis('x');
%             mx=mean(locs.xnm);
%             plot(axx,ltime,locs.xnm-mx,ltime,0*locs.time,'k', ...
%                 ltime,0*ltime+sigmax,'c',ltime,0*ltime-sigmax,'r',...
%                 ltime,0*locs.time+sxrobust,'r',ltime,0*locs.time-sxrobust,'c',...
%                 ltime,0*locs.time+sxdetrend,'m',ltime,0*locs.time-sxdetrend,'m')
%             legend(axx,'data','','std','','robust std','','detrend std','')
%             xlabel(axx,'time (ms)')
%             ylabel(axx,'x (nm)')
               title(axxy,['std(y) = ' num2str(sigmay,ff) ' nm, std(x) detrend = ' num2str(sydetrend,ff) ' nm.'])

% 
%             axy=obj.initaxis('y');
%             plot(axy,ltime,locs.ynm-mean(locs.ynm),ltime,0*locs.time,'k', ...
%                 ltime,0*locs.time+sigmay,'c',ltime,0*locs.time-sigmay,'r',...
%                 ltime,0*locs.time+syrobust,'r',ltime,0*locs.time-syrobust,'c',...
%                 ltime,0*locs.time+sydetrend,'m',ltime,0*locs.time-sydetrend,'m')
%             legend(axy,'data','','std','','robust std','','detrend std','')
%             xlabel(axy,'time (ms)')
%             ylabel(axy,'y (nm)')
%             title(axy,['std(y) = ' num2str(sigmay,ff) ' nm, std(y) robust = ' num2str(syrobust,ff) ' nm, std(y) detrend = ' num2str(sydetrend,ff) ' nm.'])

            axt=obj.setoutput('time');
            histogram(axt,dt,dtmin/2:dtmin:quantile(dt,0.995))
            hold(axt,'on')
            histogram(axt,dt,0:dtmin*0.1:quantile(dt,0.995))
            xlabel(axt,'dt (ms)')
            ylabel(axt,'frequency')
            title(axt,['dtmin = ' num2str(dtmin,ff) ' ms, dtmedian = ' num2str(dtmedian,ff) ' ms, dtmean = ' num2str(dtmean,ff) ' ms.'])
            
            axdt=obj.setoutput('dt');
            plot(axdt,ltime(2:end),dt)
            xlabel(axdt,'time (ms)')
            ylabel(axdt,'dt (ms)')

            if ~isempty(obj.locData.loc.efo)
                axe=obj.setoutput('efo');
                plot(axe,ltime,obj.locData.loc.efo(index))
                xlabel(axe,'time (ms)')
                ylabel(axe,'efo')
            end
            if ~isempty(obj.locData.loc.cfr)
                axc=obj.setoutput('cfr');
                plot(axc,ltime,obj.locData.loc.cfr(index))
                xlabel(axc,'time (ms)')
                ylabel(axc,'cfr')
            end
            if ~isempty(obj.locData.loc.eco)
                axec=obj.setoutput('eco');
                plot(axec,ltime,obj.locData.loc.eco(index))
                xlabel(axec,'time (ms)')
                ylabel(axec,'eco')
            end
            if ~isempty(obj.locData.loc.ecc)
                axcc=obj.setoutput('ecc');
                plot(axcc,ltime,obj.locData.loc.ecc(index))
                xlabel(axcc,'time (ms)')
                ylabel(axcc,'ecc')
            end

%             header=sprintf('nlocs \t on-time \t dtmin \t dtmedian \t <dt> \t sigmax \t sigmay \t sigmax robust \t sigmay robust \t sigmax detrend \t sigmay detrend \t efo med \t cfr med  \t eco med  \t ecc med  \t efc med \t filename' );
%             disp(header)
%             results=sprintf([num2str(nlocs)  '\t' num2str(ontime)  '\t' num2str(dtmin)  '\t' num2str(dtmedian) '\t' num2str(dtmean) '\t' num2str(sigmax) ...
%                  '\t' num2str(sigmay)  '\t' num2str(sxrobust)  '\t' num2str(syrobust)  '\t' num2str(sxdetrend)  '\t' num2str(sydetrend) '\t' num2str(efo) '\t' ...
%                  num2str(cfr) '\t' num2str(eco) '\t' num2str(ecc) '\t' num2str(efc) '\t' filename]  );
%             out.clipboard=results;
            out.cluster=[];

        return
        modality=p.modality.selection;
        obj.axis=obj.setoutput('profile');
        fs=obj.site.evaluation.fibrilStatistics.measurement;

        switch modality
            case 'deviation'
                dev=fs.deviation.value;
                dsmooth=fs.deviation.value_smo;
            case 'polarization'
                dev=fs.P.value;
                dsmooth=fs.P.value_smo;
        end
        if isfield(fs,'P')
            indKeptCurve = fs.P.indKeptCurve;
        else
            indKeptCurve = true(length(dev),1);
        end
        posx=(1:length(dev))'*10;
        hold(obj.axis,'off');
        plot(obj.axis,posx(indKeptCurve),dev(indKeptCurve),'-')
        hold(obj.axis,'on');
        plot(obj.axis,posx(indKeptCurve),dsmooth(indKeptCurve),'r-','LineWidth',3);
        if isfield(obj.site.evaluation,obj.name)
            out=obj.site.evaluation.(obj.name);
        else
            out.(modality).Position=[];
        end
        if isfield(obj.site.evaluation,obj.name) && isfield(obj.site.evaluation.(obj.name),modality) && ~isempty(obj.site.evaluation.(obj.name).(modality).Position)
            obj.roihandle=images.roi.Polyline(obj.axis,'Position',obj.site.evaluation.(obj.name).(modality).Position);
            out=obj.site.evaluation.(obj.name);
            addlistener(obj.roihandle,'ROIMoved',@(src,evt) obj.updateposition(src,evt));
            addlistener(obj.roihandle,'VertexAdded',@(src,evt) obj.updateposition(src,evt));
            addlistener(obj.roihandle,'VertexDeleted',@(src,evt) obj.updateposition(src,evt));
            obj.plotdistances;
        end

        end
     
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function select_callback(obj,a,b)
            if ~isempty(obj.roihandle)&&isvalid(obj.roihandle)
                delete(obj.roihandle);
            end
            obj.roihandle=drawpolyline;
%             obj.site.evaluation.(obj.name).Position=obj.roihandle.Position;
            addlistener(obj.roihandle,'ROIMoved',@(src,evt) obj.updateposition(src,evt));
            addlistener(obj.roihandle,'VertexAdded',@(src,evt) obj.updateposition(src,evt));
            addlistener(obj.roihandle,'VertexDeleted',@(src,evt) obj.updateposition(src,evt));
            updateposition(obj,a,b)
        end
        function updateposition(obj,a,b)
            if obj.getSingleGuiParameter('equaldistance')
                pos=obj.roihandle.Position;
                pos(:,1)=linspace(pos(1,1),pos(end,1),size(pos,1));
                obj.roihandle.Position=pos;
            end
            modality=obj.getSingleGuiParameter('modality').selection;
            obj.site.evaluation.(obj.name).(modality).Position=obj.roihandle.Position;
            obj.plotdistances;
        end
        function plotdistances(obj)
            modality=obj.getSingleGuiParameter('modality').selection;
            pos=obj.site.evaluation.(obj.name).(modality).Position;
            
            period2=mean(diff(pos(:,1)));
            form='%2.2f';
            period=(pos(end,1)-pos(1,1))/(size(pos,1)-1);
            
            title(obj.axis,['Period: ' num2str(period,form)])
            
            obj.site.evaluation.(obj.name).(modality).Period=period;
            
        end
%         function fit_callback(obj,a,b)
%             pos=obj.site.evaluation.(obj.name).Position;
%             period=mean(diff(pos(:,1)));
%             
%             
%         end
    end

end


function selectrange(a,b,obj)
obj.site.evaluation.(obj.name).range=obj.axstep.XLim;
obj.run(obj.getAllParameters);
end

function pard=guidef(obj)
pard.overshoott.object=struct('String','Coarsness','Style','text');
pard.overshoott.position=[1,1];
pard.overshoot.object=struct('String','1','Style','edit');
pard.overshoot.position=[1,2];
pard.overshoot.Width=0.5;


pard.fitmodet.object=struct('String','fit using','Style','text');
pard.fitmodet.position=[1,2.5];
pard.fitmode.object=struct('String',{{'mean','median'}},'Style','popupmenu');
pard.fitmode.position=[1,3.5];
pard.fitmode.Width=1.5;

pard.currentrange.object=struct('String','Current Range','Style','pushbutton','Callback',{{@selectrange,obj}});
pard.currentrange.position=[2,1];

%auto-fit


% pard.dxt.Width=3;
pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi','se_sitepixelsize'};
pard.plugininfo.type='ROI_Evaluate';
pard.plugininfo.description='';
end
