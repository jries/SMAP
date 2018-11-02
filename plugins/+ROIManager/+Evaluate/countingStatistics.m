classdef countingStatistics<interfaces.SEEvaluationProcessor
    methods
        function obj=countingStatistics(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj,p)
%             layerson=obj.locData.parameters.sr_layerson;
            out.PSFlayers=[];
            out.locplayers=[];
            out.Nlayers=[];
            roisize=obj.getPar('se_siteroi')/2;

                    [locs,ind]=obj.getLocs({'locprecnm','PSFxnm','xnm','ynm','frame'},'layer',1,'size',roisize*2,'grouping','ungrouped');  
                    [glocs,indg]=obj.getLocs({'locprecnm','PSFxnm','xnm','ynm','frame'},'layer',1,'size',roisize*2,'grouping','grouped'); 
                    out.PSF=mean(locs.PSFxnm);out.PSFg=mean(glocs.PSFxnm);
                    inc=(locs.xnm-obj.site.pos(1)).^2+(locs.ynm-obj.site.pos(2)).^2<roisize^2;
                    incg=(glocs.xnm-obj.site.pos(1)).^2+(glocs.ynm-obj.site.pos(2)).^2<roisize^2;
                    
                    
                    out.Nlocs=sum(inc);
                     out.Nlocsg=sum(incg);
                    out.locprecnm=mean(locs.locprecnm(inc));out.locprecnmg=mean(glocs.locprecnm(incg));
                    out.noblink=min(sum(diff(locs.frame(inc))>p.noblinktime)+1,out.Nlocsg);
                     out.Nlocsg=min(sum(diff(locs.frame(inc))>1)+1,out.Nlocsg);
            
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end

function pard=guidef
pard.text1.object=struct('Style','text','String','max dark time (frame)');
pard.text1.position=[2,1];
pard.text1.Width=2;
pard.noblinktime.object=struct('Style','edit','String','20');
pard.noblinktime.position=[2,3];
pard.noblinktime.Width=2;
pard.plugininfo.type='ROI_Evaluate';
end