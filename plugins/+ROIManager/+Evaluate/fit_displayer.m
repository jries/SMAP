classdef fit_displayer<interfaces.SEEvaluationProcessor

    methods
        function obj=fit_displayer(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        
        function out=run(obj,p)
            out=[];
            % set up the fitting
            m1 = functionModel('C:\Users\ries\git\ries-private\SMLMModelFitter\models\CME3DContinuous.m');         % model 1
            fitting = SMLMModelFit({m1}, 'SolverName', 'fminsearchbnd','SolverOptions',{'Display','iter'});
            % fitting = SMLMModelFit({m1}, 'SolverName', 'particleswarm','SolverOptions',{'UseVectorized',true, 'Display','iter'});
            fitting.init                     % initiation
            fitting.roiSize = 500;
            fitResult = obj.site.evaluation.SMLMModelFit.fitResult;
            fitting.allParsArg = fitResult;
            
            locs = obj.getLocs({'xnm','ynm','znm','locprecnm'},'grouping','grouped','layer',1, 'size', p.se_siteroi);
            locs.xnm = locs.xnmrot;
            locs.ynm = locs.ynmrot;
            ax1 = obj.setoutput('Projection',1);
            subplot(2,2,1,ax1);
            axes(ax1);
            ax1_2 = subplot(2,2,2);
            axes(ax1);
            ax1_3 = subplot(2,2,3);
            fitting.plot(locs,'Projection', 'xy', 'pixelSize', 2, 'axes', ax1);
            fitting.plot(locs,'Projection', 'xz', 'pixelSize', 2, 'axes', ax1_2);
            fitting.plot(locs,'Projection', 'yz', 'pixelSize', 2, 'axes', ax1_3);
            
            title(ax1, 'xy')
            title(ax1_2, 'zy')
            title(ax1_3, 'xz')
            
            ax2 = obj.setoutput('Point_cloud',1);
            xlabel(ax2, 'x')
            ylabel(ax2, 'y')
            zlabel(ax2, 'z')
            fitting.plot(locs, 'plotType', 'point', 'pixelSize', 2, 'axes',ax2);
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function pard = guidef(obj)
    pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi'};
    pard.plugininfo.type='ROI_Evaluate';
end

function filterCallback(a,b,obj)
    if obj.guihandles.filtering.Value
        obj.guihandles.lockFilter.Enable = 'on';
    else
        obj.guihandles.lockFilter.Enable = 'off';
        obj.guihandles.lockFilter.Value = 0;
    end
end