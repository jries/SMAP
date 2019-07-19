classdef VISPexport<interfaces.DialogProcessor
% Exports localizations in a format compatible with VISP: El Beheiry,
% Mohamed, and Maxime Dahan. “ViSP: Representing Single-Particle
% Localizations in Three Dimensions.” Nature Methods 10, no. 8 (August 1,
% 2013): 689–90. https://doi.org/10.1038/nmeth.2566.
    methods
        function obj=VISPexport(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;
             obj.inputParameters={'sr_layerson'};
        end
        function pard=guidef(obj)
            pard.plugininfo.type='ProcessorPlugin';
            pard.vispformat.object=struct('Style','popupmenu','String',{{'2D','2D LP','3D','3D LP'}},'Value',3);
            pard.vispformat.position=[1,1];
            pard.plugininfo.description='exports localizations in a format compatible with VISP. Exports localizations to be used with: El Beheiry, Mohamed, and Maxime Dahan. “ViSP: Representing Single-Particle Localizations in Three Dimensions.” Nature Methods 10, no. 8 (August 1, 2013): 689–90. https://doi.org/10.1038/nmeth.2566.';
        end
        
        function out=run(obj,p)
           
            [path,f]=fileparts(obj.locData.files.file(1).name);
            switch p.vispformat.selection
                case '2D'
                    ext='2d';
                case '2D LP'
                    ext='2dlp';
                case '3D'
                    ext='3d';
                case '3D LP'  
                    ext='3dlp';
            end
            [f,path]=uiputfile([path filesep f '.' ext]);
            [~,fx]=fileparts(f);
            facfwhm=2.3;
            for k=1:length(p.sr_layerson)
                if p.sr_layerson(k)
                     locs=obj.locData.getloc({'xnm','ynm','znm','locprecznm','locprecnm','phot','frame'},'position','roi','layer',k);
                     if isempty(locs.locprecznm)
                         locs.locprecznm=0*locs.xnm+1;
                     end
                     if isempty(locs.locprecnm)
                         locs.locprecnm=0*locs.xnm+1;
                     end
                    a=1;
                    switch p.vispformat.selection
                        
                        case '2D'
                            outmatrix=[locs.xnm-a*mean(locs.xnm) locs.ynm-a*mean(locs.ynm) locs.phot locs.frame];
                        case '2D LP'
                            outmatrix=[locs.xnm-a*mean(locs.xnm) locs.ynm-a*mean(locs.ynm) locs.locprecnm*facfwhm locs.locprecnm*facfwhm locs.phot locs.frame];
                        case '3D'
                            outmatrix=[locs.xnm-a*mean(locs.xnm) locs.ynm-a*mean(locs.ynm) locs.znm locs.phot locs.frame];
                        case '3D LP'  
                            outmatrix=[locs.xnm-a*mean(locs.xnm) locs.ynm-a*mean(locs.ynm) locs.znm locs.locprecnm*facfwhm locs.locprecnm*facfwhm locs.locprecznm*facfwhm   locs.phot locs.frame];
                    end
                    
                    fout=[path filesep fx '_l' num2str(k) '.' ext];
                    dlmwrite(fout,outmatrix,'delimiter','\t')
                end
            end          
            out=0;
        end
    end
end