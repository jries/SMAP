classdef clusterStatistics<interfaces.SEEvaluationProcessor
%     number of localizations, average PSF size etc.
    methods
        function obj=clusterStatistics(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj,p)
            layerson=p.sr_layerson;%obj.locData.parameters.sr_layerson;
            out.PSFlayers=[];
            out.locplayers=[];
            out.Nlayers=[];
%             roisize=obj.site.sePar.Settings.siteroi/2;
            roisize=p.se_siteroi;
            if ~p.circularroi
                roisizeh=[roisize,roisize];
            else
                roisizeh=roisize/2;
            end
            fields0={'locprecnm','xnm','ynm','znm','phot','bg','psf_nu','psf_bg'};
            fields1={'locprecznm_SALM','znm_SALM','znm_a','locprecznm_a','xnmerr','ynmerr'};
            fields=[fields0 fields1];
%             if obj.display
%              h=obj.setoutput('p1');
%              h.NextPlot='replace';
%             end
            for k=1:p.numberOfLayers
                if layerson(k)
                    locs=obj.getLocs(fields,'layer',k,'size',roisizeh);  
                    
                    for f=1:length(fields)
                        vh=locs.(fields{f});
                        vh(isnan(vh))=[];
                        vh(isinf(vh))=[];
                        
                          out.(['layers' num2str(k)]).(fields{f}).mean=mean(vh,'omitnan');
                          out.(['layers' num2str(k)]).(fields{f}).median=median(vh,'omitnan');
                          out.(['layers' num2str(k)]).(fields{f}).std=std(vh,'omitnan');
                          out.(['layers' num2str(k)]).(fields{f}).q05=quantile(vh,0.05);
                          out.(['layers' num2str(k)]).(fields{f}).q95=quantile(vh,0.95);
                          out.(['layers' num2str(k)]).(fields{f}).vals=locs.(fields{f});
                          out.(['layers' num2str(k)]).(fields{f}).NnotNAN=length(vh);
                    end
                    
                     out.(['layers' num2str(k)]).Nlocs=length(locs.xnm);
                     
                end
            end
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end

function pard=guidef
pard.circularroi.object=struct('Style','checkbox','String','Use circular ROI');
pard.circularroi.position=[1,1];
pard.circularroi.Width=4;
pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi'};
pard.plugininfo.type='ROI_Evaluate';
pard.plugininfo.description='number of localizations, average PSF size etc.';
end