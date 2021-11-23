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