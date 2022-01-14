classdef GetExposureTime<interfaces.DialogProcessor&interfaces.SEProcessor

    properties
    end
    methods
        function obj=GetExposureTime(varargin)        
            obj@interfaces.DialogProcessor(varargin{:});
        end
        
        function out=run(obj,p)  
            exposuretime=obj.locData.files.file.info.exposure;
            filename=obj.locData.files.file.name;
            out.exposuretime=exposuretime;
            out.filename=filename;
            
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end

    end
end






function pard=guidef(obj)
pard.t1.object=struct('String','Get exposure time','Style','text');
pard.t1.position=[1,1];
% 
% pard.corners.object=struct('String','8','Style','edit');
% pard.corners.position=[1,2];
% pard.corners.Width=0.5;
% 
% pard.t2.object=struct('String','Proteins/Corner','Style','text');
% pard.t2.position=[1,3];
% 
% pard.rings.object=struct('String','4','Style','edit');
% pard.rings.position=[1,4];
% pard.rings.Width=0.5;
% 
% 
% pard.psfcheck.object=struct('String','PSF range','Style','checkbox');
% pard.psfcheck.position=[2,1];
% 
% pard.PSFrange.object=struct('String','80 150','Style','edit');
% pard.PSFrange.position=[2,2];
% pard.PSFrange.Width=1;
% 
% pard.t3.object=struct('String','fit range histogram','Style','text');
% pard.t3.position=[2,3];
% 
% pard.fitrange.object=struct('String','3 8','Style','edit');
% pard.fitrange.position=[2,4];
% pard.fitrange.Width=1;
% 
% pard.filecheck.object=struct('String','filenumbers','Style','checkbox');
% pard.filecheck.position=[3,1];
% 
% pard.filenumbers.object=struct('String','1:100','Style','edit');
% pard.filenumbers.position=[3,2];
% pard.filenumbers.Width=1;
% 
% % pard.bootstrap.object=struct('String','Bootstrap error bars for time analysis','Style','checkbox');
% % pard.bootstrap.position=[4,1];
% % pard.bootstrap.Width=2;
% 
% pard.copy2page.object=struct('String','Copy to own page','Style','checkbox');
% pard.copy2page.position=[6,1];
% pard.copy2page.Width=2;

pard.plugininfo.type='ROI_Analyze';
% pard.plugininfo.description='Calculates the absolute effective labeling efficiency from NPC images using the results of evaluator: NPCLabelingQuantify_s. See:Thevathasan, Jervis Vermal, Maurice Kahnwald, Konstanty Cieśliński, Philipp Hoess, Sudheer Kumar Peneti, Manuel Reitberger, Daniel Heid, et al. “Nuclear Pores as Versatile Reference Standards for Quantitative Superresolution Microscopy.�? BioRxiv, March 20, 2019, 582668. https://doi.org/10.1101/582668.';

end