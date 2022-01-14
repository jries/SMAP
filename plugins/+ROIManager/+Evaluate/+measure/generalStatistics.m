classdef generalStatistics<interfaces.SEEvaluationProcessor
%     number of localizations, average PSF size etc.
    methods
        function obj=generalStatistics(varargin)        
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
%             if obj.display
%              h=obj.setoutput('p1');
%              h.NextPlot='replace';
%             end
            for k=1:p.numberOfLayers
                if layerson(k)
%                     locs=obj.getLocs({'locprecnm','PSFxnm','xnm','ynm','phot','bg','numberInGroup'},'layer',k,'size','freeroi');  
                    locs=obj.getLocs({'locprecnm','PSFxnm','xnm','ynm','phot','bg','numberInGroup'},'layer',k,'size',roisizeh);  
%                    figure(88);plot(locs.xnm,locs.ynm,'.')
                    if isfield(locs,'PSFxnm')
                    psf=mean(locs.PSFxnm);
                    out.PSFlayers(end+1)=psf;
                    out.(['layers' num2str(k)]).PSF=psf;
                    end
                    
                    locp=median(locs.locprecnm);
                    out.locplayers(end+1)=locp;
                    out.(['layers' num2str(k)]).locp=locp;
                    
                    N=length(locs.xnm);
                    out.Nlayers(end+1)=N;
                    out.(['layers' num2str(k)]).N=N;
                    
                    out.(['layers' num2str(k)]).meanx=mean(locs.xnm);
                    out.(['layers' num2str(k)]).meany=mean(locs.ynm); 
                    out.(['layers' num2str(k)]).medianx=median(locs.xnm);
                    out.(['layers' num2str(k)]).mediany=median(locs.ynm);
                    out.(['layers' num2str(k)]).meanphot=mean(locs.phot);
                    out.(['layers' num2str(k)]).medianphot=median(locs.phot);
                    out.(['layers' num2str(k)]).meanbg=mean(locs.bg);
                    out.(['layers' num2str(k)]).medianbg=median(locs.bg);
                    out.(['layers' num2str(k)]).meanlifetime=mean(locs.numberInGroup);
                    if ~isempty(locs.bg)
                    out.(['layers' num2str(k)]).maxbg=mode(round(locs.bg));
                    else
                    out.(['layers' num2str(k)]).maxbg=NaN;
                    end
                end
            end
            
            out.Nch=[];
            out.Nchg=[];
            
            maxc=4;
            locs=obj.getLocs({'channel'},'size',roisizeh,'grouping','ungrouped'); 
            glocs=obj.getLocs({'channel'},'size',roisizeh,'grouping','grouped');           
            for c=0:maxc
                ind=locs.channel==c;
                Nch=sum(ind);
                if Nch>0
                    out.Nch(end+1)=Nch;
                    out.(['channels' num2str(c)]).N=Nch;
                    
                    Nchg=sum(glocs.channel==c);
                    out.Nchg(end+1)=Nchg;
                    out.(['channels' num2str(c)]).Ng=Nchg;
                    
                   
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