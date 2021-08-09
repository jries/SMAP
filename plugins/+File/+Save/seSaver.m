classdef seSaver<interfaces.DialogProcessor
    % this is for alpha-syn project
    methods
        function obj=seSaver(varargin)  
                obj@interfaces.DialogProcessor(varargin{:}) ;
                obj.inputParameters={'filelist_long','mainfile','mainGui','numberOfLayers','sr_layerson','sr_pixrec','layers','sr_image','sr_pos','group_dt','group_dx'};
        end
        
        function out=save(obj,p)
            obj.status('save SE file')
            lastSMLFile = obj.getPar('lastSMLFile');
            defaultFn = replace(lastSMLFile, '_sml', '_se');
            
            % only used sites will be saved
            lUse = getFieldAsVector(obj.locData.SE.sites,'annotation.use');
            subSites = obj.locData.SE.sites(lUse);
            fibrilStatistics = getFieldAsVector(subSites,'evaluation.fibrilStatistics');
            fibrilDynamics = getFieldAsVector(subSites,'evaluation.fibrilDynamics');
            fibrilStraightener = getFieldAsVector(subSites,'evaluation.fibrilStraightener');
            annotatePeaks = getFieldAsVector(subSites,'evaluation.annotatePeaks');
            fnMeasurement = {'deviation','P','intensity','intensity_rescaled','P_m','deviation_m'};
            for k = 1:length(fibrilStatistics)
                singleSites{k}.pos = subSites(k).pos;
                singleSites{k}.annotation = subSites(k).annotation;
                singleSites{k}.fibrilStatistics = fibrilStatistics{k};
                singleSites{k}.fibrilStatistics.kymograph = [];
                singleSites{k}.fibrilStatistics.GuiParameters = [];
                for l = 1:length(fnMeasurement)
                    singleSites{k}.fibrilStatistics.measurement.(fnMeasurement{l}).raw = [];
                    singleSites{k}.fibrilStatistics.measurement.(fnMeasurement{l}).fft = [];
                end
                singleSites{k}.fibrilDynamics = fibrilDynamics{k};
                singleSites{k}.fibrilDynamics.GuiParameters = [];
                singleSites{k}.fibrilStraightener = fibrilStraightener{k};
                singleSites{k}.fibrilStraightener.indFibrilLocs = [];
                singleSites{k}.annotatePeaks = annotatePeaks{k};
            end
            uisave('singleSites', defaultFn)
            obj.status('save done')
            out = [];
        end
        function pard=guidef(obj)
           pard.plugininfo.type='SaverPlugin';
        end
        function out = run(obj,p)
            obj.save(p)
            out = [];
        end        

    end
end
