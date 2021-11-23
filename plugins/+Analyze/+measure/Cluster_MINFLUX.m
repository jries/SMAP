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
            locs=obj.locData.getloc({'xnm','ynm','time','numberInGroup','filenumber'},'layer',find(obj.getPar('sr_layerson')),'Position','roi');
            filelist=obj.getPar('filelist_short');
            filename=filelist.String{mode(locs.filenumber)};
            dt=diff(locs.time);
            dtmin=min(dt);
            dtmedian=median(dt);
            dtmean=mean(dt);

            sigmax=std(locs.xnm);sigmay=std(locs.ynm);
            sxdetrend=std(diff(locs.xnm))/sqrt(2);sydetrend=std(diff(locs.ynm))/sqrt(2);
            [~, sxrobust]=robustMean(locs.xnm); [~, syrobust]=robustMean(locs.ynm);
            
            %graphs
                 ff='%2.1f';
            axx=obj.initaxis('x');
            mx=mean(locs.xnm);
            plot(axx,locs.time,locs.xnm-mx,locs.time,0*locs.time,'k', ...
                locs.time,0*locs.time+sigmax,'c',locs.time,0*locs.time-sigmax,'r',...
                locs.time,0*locs.time+sxrobust,'r',locs.time,0*locs.time-sxrobust,'c',...
                locs.time,0*locs.time+sxdetrend,'m',locs.time,0*locs.time-sxdetrend,'m')
            legend(axx,'data','','std','','robust std','','detrend std','')
            xlabel(axx,'time (ms)')
            ylabel(axx,'x (nm)')
               title(axy,['std(x) = ' num2str(sigmax,ff) ' nm, std(x) robust = ' num2str(sxrobust,ff) ' nm, std(x) detrend = ' num2str(sxdetrend,ff) ' nm.'])


            axy=obj.initaxis('y');
            plot(axy,locs.time,locs.ynm-mean(locs.ynm),locs.time,0*locs.time,'k', ...
                locs.time,0*locs.time+sigmay,'c',locs.time,0*locs.time-sigmay,'r',...
                locs.time,0*locs.time+syrobust,'r',locs.time,0*locs.time-syrobust,'c',...
                locs.time,0*locs.time+sydetrend,'m',locs.time,0*locs.time-sydetrend,'m')
            legend(axy,'data','','std','','robust std','','detrend std','')
            xlabel(axy,'time (ms)')
            ylabel(axy,'y (nm)')
            title(axy,['std(y) = ' num2str(sigmay,ff) ' nm, std(y) robust = ' num2str(syrobust,ff) ' nm, std(y) detrend = ' num2str(sydetrend,ff) ' nm.'])

            axt=obj.initaxis('time');
            histogram(axt,dt,dtmin/2:dtmin:quantile(dt,0.995))
            hold(axt,'on')
            histogram(axt,dt,0:dtmin*0.1:quantile(dt,0.995))
            xlabel(axt,'dt (ms)')
            ylabel(axt,'frequency')
       
            title(axt,['dtmin = ' num2str(dtmin,ff) ' ms, dtmedian = ' num2str(dtmedian,ff) ' ms, dtmean = ' num2str(dtmean,ff) ' ms.'])
 
            header=sprintf('filename \t dtmin \t dtmedian \t dtmean \t sigmax \t sigmay \t sigmaxrobust \t sigmayrobust \t sigmaxdetrend \t sigmaydetrend' );
            disp(header)
            results=sprintf([filename '\t' num2str(dtmin)  '\t' num2str(dtmedian) '\t' num2str(dtmean) '\t' num2str(sigmax) ...
                 '\t' num2str(sigmay)  '\t' num2str(sxrobust)  '\t' num2str(syrobust)  '\t' num2str(sxdetrend)  '\t' num2str(sydetrend)]  );
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