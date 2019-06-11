classdef UseSitesInROI<interfaces.DialogProcessor&interfaces.SEProcessor
    methods
        function obj=UseSitesInROI(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'se_sitefov','se_cellpixelsize','se_siteroi'};
        end
        
        function out=run(obj,p)  
          hroi=obj.getPar('sr_roihandle');
          imbw=hroi.createMask;
          pos=getFieldAsVectorInd(obj.SE.sites,'pos');
          pr=obj.getPar('sr_pixrec');
          srpos=obj.getPar('sr_pos');
          srsize=obj.getPar('sr_size');
          prel=(pos(:,1:2)-(srpos-srsize))/pr;
          ing=withinmask(imbw,prel(:,1),prel(:,2));
          for k=1:length(obj.SE.sites)
              obj.SE.sites(k).annotation.use=ing(k);
          end
          out=[];
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end




function pard=guidef

pard.plugininfo.type='ROI_Analyze';
end
