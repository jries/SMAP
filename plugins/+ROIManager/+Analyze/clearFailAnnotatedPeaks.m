classdef clearFailAnnotatedPeaks<interfaces.DialogProcessor&interfaces.SEProcessor
    methods
        function obj=clearFailAnnotatedPeaks(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'se_viewer', 'se_siteroi', 'se_sitefov'};
            obj.showresults=true;
        end
        
        function out=run(obj,p)
            se = obj.locData.SE;
            sites = se.sites;
            list1 = getFieldAsVector(sites,'annotation.list1.value');
            for k = 1:se.numberOfSites
                annotatePeaks = sites(k).evaluation.annotatePeaks;
                switch list1(k)
                    case 3
                        if isfield(annotatePeaks, 'deviation')
                            sites(k).evaluation.annotatePeaks = rmfield(annotatePeaks, 'deviation');
                        end
                    case 4
                        if isfield(annotatePeaks, 'polarization')
                            sites(k).evaluation.annotatePeaks = rmfield(annotatePeaks, 'polarization');
                        end
                    case 5
                        if isfield(annotatePeaks, 'deviation')
                            sites(k).evaluation.annotatePeaks = rmfield(annotatePeaks, 'deviation');
                        end
                        if isfield(annotatePeaks, 'polarization')
                            sites(k).evaluation.annotatePeaks = rmfield(annotatePeaks, 'polarization');
                        end
                    otherwise
                end
            end
            out = [];
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function pard=guidef(obj)

pard.plugininfo.type='ROI_Analyze';

end