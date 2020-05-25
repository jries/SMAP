classdef countingStatistics_kinetochore<interfaces.SEEvaluationProcessor
%     calcuculate number of localizations with different dark times. Used
%     for locsfromSE plugin (see counting)
    methods
        function obj=countingStatistics_kinetochore(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj,p)
%             layerson=obj.locData.parameters.sr_layerson;
            out.PSFlayers=[];
            out.locplayers=[];
            out.Nlayers=[];
            roisize=obj.getPar('se_siteroi')/2;

                    [locs,ind]=obj.getLocs({'locprecnm','PSFxnm','xnm','ynm','frame'},'layer',1,'size','freeroi','grouping','ungrouped');  
                    [glocs,indg]=obj.getLocs({'locprecnm','PSFxnm','xnm','ynm','frame'},'layer',1,'size','freeroi','grouping','grouped'); 
                    out.PSF=mean(locs.PSFxnm);out.PSFg=mean(glocs.PSFxnm);
                    inc=(locs.xnm-obj.site.pos(1)).^2+(locs.ynm-obj.site.pos(2)).^2<roisize^2;
                    incg=(glocs.xnm-obj.site.pos(1)).^2+(glocs.ynm-obj.site.pos(2)).^2<roisize^2;
                    % need to work on this more
                    structureType = obj.site.annotation.list4.string{obj.site.annotation.list4.value};
                    copynumber = obj.site.annotation.list3.string{obj.site.annotation.list3.value};
                    
                    if isempty(obj.site.annotation.line3)
                        areaMask = roisize^2;
                    else
                        xMask = obj.site.annotation.line3.pos(:,1);
                        yMask = obj.site.annotation.line3.pos(:,2);
                        areaMask = polyarea(xMask, yMask)*1e6;
                    end
                    
                    if strcmp(structureType, 'kinetochore')&&strcmp(copynumber, '32 copies')
                        out.Nlocs=sum(inc)/2;
                        out.Nlocsg=sum(incg)/2;
                        out.areaMask=areaMask/2;
                    else
                        out.Nlocs=sum(inc);
                        out.Nlocsg=sum(incg);
                        out.areaMask=areaMask;
                    end
                    
                    
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
pard.plugininfo.description='calcuculate number of localizations with different dark times. Used for locsfromSE plugin (see counting)';
end