classdef Cluster_MINFLUX<interfaces.DialogProcessor
    % LINEPROFILE Calculates profiles along a linear ROI and fits it with a
    % model of choice. Flat: step function convolved with Gaussian
    % (=Erf). Disk: Projection of a homogeneously filled disk, convolved
    % with Gaussian. Ring: Projection of a ring, convolved with
    % Gaussian. Distance: Two Gaussians in a distance d.
    methods
        function obj=Cluster_MINFLUX(varargin)        
            obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'sr_layerson','linewidth_roi','znm_min','znm_max','sr_pixrec','layernames'};
            obj.showresults=true;
        end
        
        function out=run(obj,p)      
            locs=obj.locData.getloc({'xnm','ynm','znm','time','numberInGroup','filenumber','efo','cfr','eco','ecc','efc','id'},'layer',find(obj.getPar('sr_layerson')),'Position','roi');
            filelist=obj.getPar('filelist_short');
            filename=filelist.String{mode(locs.filenumber)};
            dt=diff(locs.time);
            dtmin=min(dt(dt>0));
            dtmedian=median(dt);
            dtmean=mean(dt);
            efo=median(locs.efo,'omitnan');
            cfr=median(locs.cfr,'omitnan');
            eco=median(locs.eco,'omitnan');
            ecc=median(locs.ecc,'omitnan');
            efc=median(locs.efc,'omitnan');
            nlocs=length(locs.time);
            ontime=max(locs.time)-min(locs.time);

            ltime=locs.time-min(locs.time);

            sigmax=std(locs.xnm);sigmay=std(locs.ynm);
            sxdetrend=std(diff(locs.xnm))/sqrt(2);sydetrend=std(diff(locs.ynm))/sqrt(2);
            [~, sxrobust]=robustMean(locs.xnm); [~, syrobust]=robustMean(locs.ynm);
            
            %graphs
                 ff='%2.1f';
            axx=obj.initaxis('x');
            mx=mean(locs.xnm);
            plot(axx,ltime,locs.xnm-mx,ltime,0*locs.time,'k', ...
                ltime,0*ltime+sigmax,'c',ltime,0*ltime-sigmax,'r',...
                ltime,0*locs.time+sxrobust,'r',ltime,0*locs.time-sxrobust,'c',...
                ltime,0*locs.time+sxdetrend,'m',ltime,0*locs.time-sxdetrend,'m')
            legend(axx,'data','','std','','robust std','','detrend std','')
            xlabel(axx,'time (ms)')
            ylabel(axx,'x (nm)')
               title(axx,['std(x) = ' num2str(sigmax,ff) ' nm, std(x) robust = ' num2str(sxrobust,ff) ' nm, std(x) detrend = ' num2str(sxdetrend,ff) ' nm.'])


            axy=obj.initaxis('y');
            plot(axy,ltime,locs.ynm-mean(locs.ynm),ltime,0*locs.time,'k', ...
                ltime,0*locs.time+sigmay,'c',ltime,0*locs.time-sigmay,'r',...
                ltime,0*locs.time+syrobust,'r',ltime,0*locs.time-syrobust,'c',...
                ltime,0*locs.time+sydetrend,'m',ltime,0*locs.time-sydetrend,'m')
            legend(axy,'data','','std','','robust std','','detrend std','')
            xlabel(axy,'time (ms)')
            ylabel(axy,'y (nm)')
            title(axy,['std(y) = ' num2str(sigmay,ff) ' nm, std(y) robust = ' num2str(syrobust,ff) ' nm, std(y) detrend = ' num2str(sydetrend,ff) ' nm.'])
            
            if ~isempty(locs.znm)
                sigmaz=std(locs.znm);
                szdetrend=std(diff(locs.znm))/sqrt(2);
                [~, szrobust]=robustMean(locs.znm); 
                axz=obj.initaxis('z');
                plot(axz,ltime,locs.znm-mean(locs.znm),ltime,0*locs.time,'k', ...
                    ltime,0*locs.time+sigmaz,'c',ltime,0*locs.time-sigmaz,'r',...
                    ltime,0*locs.time+szrobust,'r',ltime,0*locs.time-szrobust,'c',...
                    ltime,0*locs.time+szdetrend,'m',ltime,0*locs.time-szdetrend,'m')
                legend(axz,'data','','std','','robust std','','detrend std','')
                xlabel(axz,'time (ms)')
                ylabel(axz,'z (nm)')
                title(axz,['std(z) = ' num2str(sigmaz,ff) ' nm, std(z) robust = ' num2str(szrobust,ff) ' nm, std(z) detrend = ' num2str(szdetrend,ff) ' nm.'])
                tz='nlocs \t on-time \t dtmin \t dtmedian \t <dt> \t sigmax \t sigmay \t sigmaz \t sigmax robust \t sigmay robust \t sigmaz robust \t sigmax detrend \t sigmay detrend \t sigmaz detrend \t efo med \t cfr med  \t eco med  \t ecc med  \t efc med \t filename';
                tsig=['\t' num2str(sigmax) '\t' num2str(sigmay) '\t' num2str(sigmaz)  '\t' num2str(sxrobust)  '\t' num2str(syrobust) '\t' num2str(szrobust)  '\t' num2str(sxdetrend)  '\t' num2str(sydetrend) '\t' num2str(szdetrend)];
            else
                tz='nlocs \t on-time \t dtmin \t dtmedian \t <dt> \t sigmax \t sigmay \t sigmax robust \t sigmay robust \t sigmax detrend \t sigmay detrend \t efo med \t cfr med  \t eco med  \t ecc med  \t efc med \t filename' ;
                tsig=['\t' num2str(sigmax) '\t' num2str(sigmay)  '\t' num2str(sxrobust)  '\t' num2str(syrobust)  '\t' num2str(sxdetrend)  '\t' num2str(sydetrend)];
            end

            axt=obj.initaxis('time');
            histogram(axt,dt,dtmin/2:dtmin:max(quantile(dt,0.995),dtmin/2+2*dtmin))
            hold(axt,'on')
            histogram(axt,dt,0:dtmin*0.1:quantile(dt,0.995))
            xlabel(axt,'dt (ms)')
            ylabel(axt,'frequency')
            title(axt,['dtmin = ' num2str(dtmin,ff) ' ms, dtmedian = ' num2str(dtmedian,ff) ' ms, dtmean = ' num2str(dtmean,ff) ' ms.'])
            
            axdt=obj.initaxis('dt');
            plot(axdt,ltime(2:end),dt)
            xlabel(axdt,'time (ms)')
            ylabel(axdt,'dt (ms)')

            if ~isempty(locs.efo)
                axe=obj.initaxis('efo');
                plot(axe,ltime,locs.efo)
                xlabel(axe,'time (ms)')
                ylabel(axe,'efo')
            end
            if ~isempty(locs.cfr)
                axc=obj.initaxis('cfr');
                plot(axc,ltime,locs.cfr)
                xlabel(axc,'time (ms)')
                ylabel(axc,'cfr')
            end
            if ~isempty(locs.eco)
                axec=obj.initaxis('eco');
                plot(axec,ltime,locs.eco)
                xlabel(axec,'time (ms)')
                ylabel(axec,'eco')
            end
            if ~isempty(locs.ecc)
                axcc=obj.initaxis('ecc');
                plot(axcc,ltime,locs.ecc)
                xlabel(axcc,'time (ms)')
                ylabel(axcc,'ecc')
            end

            header=sprintf(tz);
            disp(header)
            results=sprintf([num2str(nlocs)  '\t' num2str(ontime)  '\t' num2str(dtmin)  '\t' num2str(dtmedian) '\t' num2str(dtmean) tsig '\t' num2str(efo) '\t' ...
                 num2str(cfr) '\t' num2str(eco) '\t' num2str(ecc) '\t' num2str(efc) '\t' filename]  );
            out.clipboard=results;

        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end




function pard=guidef(obj)

%  p(1).value=0; p(1).on={}; p(1).off={'binwidth'};
% p(2).value=1; p(2).on={'binwidth'}; p(2).off={};
% pard.setbinwidth.object=struct('String','set binwidth (nm) (otherwise: pixelsize):','Style','checkbox','Value',1,'Callback',{{@obj.switchvisible,p}});
% pard.setbinwidth.position=[1,1];
% pard.setbinwidth.Width=3;
% 
% pard.binwidth.object=struct('String','2','Style','edit');
% pard.binwidth.position=[1,3.5];
% pard.binwidth.Width=0.5;
% pard.binwidth.TooltipString='Binwidth for profiles. If not checked, use pixel size of reconstruction';
% pard.setbinwidth.TooltipString=pard.binwidth.TooltipString;



pard.plugininfo.description=sprintf('');
pard.plugininfo.type='ProcessorPlugin';
end