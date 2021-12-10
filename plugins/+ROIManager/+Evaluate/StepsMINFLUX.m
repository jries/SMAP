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
            

            if p.splitmerge
                stairs(ax2,tv(inds),mv,'m')
                if isempty(p.splitmergestep)
                    stepsize=median(steps.properties.StepSize);
                else 
                    stepsize=p.splitmergestep;
                end
                [indstep,stepvalue,posy]=splitmergefit(xr(indtime),stepsize,p,steps,yr(indtime));
                %stairs(ax2,tv(indstep),stepvalue,'r')
                inds=indstep;
                mv=stepvalue;
            end
            stairs(ax2,[tv(inds) ;tv(end)],[mv ;mv(end)],'r')


            out.steps.indstep=indstep;
            out.steps.steptime=tv(indstep);
            out.steps.stepvalue=stepvalue;
            out.steps.stepsize=diff(stepvalue);
            out.steps.possteps.x=stepvalue;
            out.steps.possteps.y=posy;
            out.steps.possteps.time=tv(indstep);

          goff=median(mod(stepvalue,16),'omitnan');
          axxy=obj.setoutput('xy');
            hold(axxy,'off')
           plot(axxy, xr-goff, yr)
           axis(axxy,'equal')
           xlabel(axxy,'x (nm)')
           ylabel(axxy,'y (nm)')
            hold(axxy,'on')
            scatter(axxy,stepvalue-goff,posy,'k')
            grid(axxy,'on')
            axm=-16:-16:axxy.XLim(1);
            axxy.XTick=[axm(end:-1:1) 0:16:axxy.XLim(2)];
            axxy.YTick=round((axxy.YLim(1):6:axxy.YLim(2))/6)*6;

%             axxy.XTick=goff:16:max(axxy.XTick);
            
%             out.steps.meanpos  calculate mean postions of step in x, y

            

            dmv=diff(mv);
            for k=1:length(dmv)
                text(ax2,tv(indstep(k+1)),mean(mv(k:k+1)),num2str(dmv(k),'%2.0f'),'FontSize',11,'Color','magenta')
            end
            obj.axstep=ax2;

            ax3=obj.setoutput('stephist');
            %stepsize=steps.properties.StepSize;
            stepsize=dmv;
            histogram(ax3,stepsize,min(stepsize):.5:max(stepsize));
            out.stepfinder=steps;

            %step time
            steptime=diff(tv(inds));
            out.steps.dwelltime=steptime;
            axstept=obj.setoutput('steptime');
            dt=max(0.1,round(mean(steptime)/5));
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
%             out.cluster=[];
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end     
    end

end


function selectrange(a,b,obj)
obj.site.evaluation.(obj.name).range=obj.axstep.XLim;
obj.run(obj.getAllParameters);
end

function [istepfit,svalfit,svalfity]=splitmergefit(x,stepsize,p,steps,y)
if contains(p.fitmode.selection,'mean')
    mfun=@mean;
else
    mfun=@median;
end

istep=[1 ; steps.properties.IndexStep];
sval=[steps.properties.LevelBefore; steps.properties.LevelAfter(end)];
    svalfit=sval;
    istepfit=istep;
for s=1:3
    [istepfit, svalfit]=mergesplit(istepfit,svalfit,stepsize);
    % recalculate sval based on x and istep2

    for k=1:10
        istepfit=moveind(x,svalfit,istepfit);
        svalfit=stepvalue(x,istepfit,mfun);
    end
    
end

if nargin>3 %y passed on
    svalfity=stepvalue(y,istepfit,mfun);
else 
    svalfity=[];
end
end

function sval=stepvalue(x,istep,fun)
istep=[istep ;length(x)+1];
for k=length(istep)-1:-1:1
    sval(k,1)=fun(x(istep(k):istep(k+1)-1));
end
end

function istep=moveind(x,sval,istep)
%step=[istep; length(x)];
d=3;
for k=2:length(istep)-1
    ih=istep(k);
    if ih>d+1
        ind=find(abs(x(ih-d:ih-1)-sval(k))<abs(x(ih-d:ih-1)-sval(k-1)));
    else
        ind=[];
    end
    if ~isempty(ind)
        istep(k)=istep(k)-d-1+min(ind);
    else
    ind=find(abs(x(ih+1:ih+d)-sval(k))>abs(x(ih+1:ih+d)-sval(k-1)));
    if ~isempty(ind)
       istep(k)=istep(k)+max(ind);
    end
    end

%      if abs(x(ih-1)-sval(k))<abs(x(ih-1)-sval(k-1))
%          istep(k)=istep(k)-1;
%      elseif abs(x(ih+1)-sval(k))<abs(x(ih+1)-sval(k+1))
%          istep(k)=istep(k)+1;
%      end
end
%istep=istep(1:end-1);
end

function [istep, sval]=mergesplit(istep,sval,stepsize)
%figure(88);hold off;
%stairs(istep,sval); 
stepv=diff(sval);
step2=stepv(1:end-1)+stepv(2:end);
sstep=find(abs(step2)<stepsize*1.4 & abs(stepv(1:end-1))<stepsize*.7);
k=1;
sstepc=sstep;
while k<length(sstepc)
    ind=find(sstepc==sstepc(k)+1);
    if ~isempty(ind)
        sstepc(ind)=[];
        %k=k+1;
    end
    k=k+1;
end
%sstep=abs(diff(sval))<stepsize*0.6;
%smallstep=find(diff(sstep).*sstep(1:end-1));
%indout=(smallstep);
istep(sstepc+2)=round((istep(sstepc+1)+istep(sstepc+2))/2);
istep(sstepc+1)=[];
%sval(sstep+2)=sval(sstep+3);
sval(sstepc+1)=[];

%stairs(istep,sval);

bigstep=find(abs(diff(sval))>stepsize*1.4)+1; %later pass on min / max step size

for k=1:length(bigstep)
    istep=[istep(1:bigstep(k)) ;istep(bigstep(k)) ;istep(bigstep(k)+1:end)];
    sval=[sval(1:bigstep(k)-1) ;(sval(bigstep(k)-1)+sval(bigstep(k)))/2 ;sval(bigstep(k):end)];
    
    istep(bigstep(k))=istep(bigstep(k))-1;
    istep(bigstep(k)+1)=istep(bigstep(k)+1)+1;
    bigstep=bigstep+1;
end
%stairs(istep,sval);
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
pard.currentrange.Width=2;
%auto-fit
p(1).value=0; p(1).on={}; p(1).off={'splitmerget','splitmergestep'};
p(2).value=1; p(2).on=p(1).off; p(2).off={};
pard.splitmerge.object=struct('Value',1,'String','Split and merge','Style','checkbox','Callback',{{@obj.switchvisible,p}});
pard.splitmerge.position=[3,1];
pard.splitmerge.Width=2;
pard.splitmerget.object=struct('String','step','Style','text');
pard.splitmerget.position=[3,3];
pard.splitmergestep.object=struct('String','','Style','edit');
pard.splitmergestep.position=[3,4];

% pard.dxt.Width=3;
pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi','se_sitepixelsize'};
pard.plugininfo.type='ROI_Evaluate';
pard.plugininfo.description='';
end
