classdef CME2DfittingUpdate<interfaces.SEEvaluationProcessor
    methods
        function obj=CME2DfittingUpdate(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        
        function out=run(obj,p)
            p.setMargin = 30;
            roisize = obj.getPar('se_siteroi');
            out = [];
            
            if ~isfield(obj.site.evaluation, 'CME2DfittingUpdate')
                out.modified = 0;
                obj.site.evaluation.CME2DfittingUpdate = out;
                
            end
            if ~obj.site.evaluation.CME2DfittingUpdate.modified
                % update position

                refyRoi = roisize/2-p.setMargin;
                refxRoi = 0;

                refxSite = obj.site.evaluation.CME2Dfitting2.parameters(1)-500/2;
                refySite = obj.site.evaluation.CME2Dfitting2.parameters(2)-500/2;
                syRoiRot = refxSite-refyRoi;
                sxRoiRot = refySite-refxRoi;
                [sxRoi, syRoi] = rotcoord(sxRoiRot, syRoiRot, (-(obj.site.annotation.rotationpos.angle)/180)*pi);
                obj.site.pos = obj.site.pos + [sxRoi syRoi 0];

                % update angle
                siteRotPos = obj.site.annotation.rotationpos.pos;
                [obj.site.annotation.rotationpos.pos(:,1),obj.site.annotation.rotationpos.pos(:,2)] = rotcoord(siteRotPos(:,1), siteRotPos(:,2), obj.site.evaluation.CME2Dfitting2.parameters(6)*pi/180);
                out.modified = 1;
             else
                disp('This site had been modified by CME2DfittingUpdate already.')
                out.modified = 1;
             end

        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end




function pard=guidef

pard.t1.object=struct('String','average sites','Style','text');
pard.t1.position=[1,1];
pard.t1.Width=1;

pard.post.object=struct('String','all or annotated:use','Style','text');
pard.post.position=[2,1];
pard.post.Width=1;

% pard.pos.object=struct('String','0,0','Style','edit');
% pard.pos.position=[2,2];
% pard.pos.Width=1;


pard.sortselection.object=struct('String',{{'all','use'}},'Style','popupmenu');
pard.sortselection.position=[2,2];
pard.sortselection.Width=1;

pard.namet.object=struct('String','name','Style','text');
pard.namet.position=[3,1];
pard.namet.Width=1;

pard.name.object=struct('String','average','Style','edit');
pard.name.position=[3,2];
pard.name.Width=1;

pard.addfile.object=struct('String','add average as new data set','Style','checkbox', 'Value', 1);
pard.addfile.position=[3,3];
pard.addfile.Width=2;

pard.t4.object=struct('String','Range','Style','text');
pard.t4.position=[4,1];
pard.t4.Width=1;

pard.siteRange.object=struct('String','Start binSize end','Style','edit');
pard.siteRange.position=[4,2];
pard.siteRange.Width=1;

pard.t5.object=struct('String','Row of sites','Style','text');
pard.t5.position=[5,1];
pard.t5.Width=1;

pard.rowSites.object=struct('String','1','Style','edit');
pard.rowSites.position=[5,2];
pard.rowSites.Width=1;

pard.inOneFile.object=struct('String','same file','Style','checkbox', 'Value', 1);
pard.inOneFile.position=[5,3];
pard.inOneFile.Width=2;

pard.plugininfo.type='ROI_Evaluate';


end
