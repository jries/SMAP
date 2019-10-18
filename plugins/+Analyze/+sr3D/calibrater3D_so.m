classdef calibrater3D_so<interfaces.DialogProcessor
%     Starts the gui to calibrate an experimental PSF model from bead stacks
% from: Li, Yiming, Markus Mund, Philipp Hoess, Joran Deschamps, Ulf Matti,
% Bianca Nijmeijer, Vilma Jimenez Sabinina, Jan Ellenberg, Ingmar Schoen,
% and Jonas Ries. “Real-Time 3D Single-Molecule Localization Using
% Experimental Point Spread Functions.” Nature Methods 15, no. 5 (April 9,
% 2018): 367–69. https://doi.org/10.1038/nmeth.4661.
    properties
        SXY
    end
    methods
        function obj=calibrater3D_so(varargin)   
            obj@interfaces.DialogProcessor(varargin{:}) ;
         end
        
        function initGui(obj)
        end
        function out=run(obj,p)
            
%             fit3ddir=strrep(pwd,'SMAP','fit3D');
%             if exist(fit3ddir,'file') && ~isdeployed
%                 addpath(fit3ddir);
%             end
%             
%             fit3ddir=strrep(pwd,'SMAP','fit3Dcspline');
%             if exist(fit3ddir,'file') && ~isdeployed
%                 addpath(fit3ddir);
%             end          
            fit3ddir='fit3Dcspline';
            smapdir=pwd;
             smappos.P=obj.P;

             smappos.fit3ddir=fit3ddir;
             smappos.smapdir=smapdir;
             cg=calibrate3D_GUI_g(smappos);    
             out=[];
        end
        
 
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end


function pard=guidef(obj)
pard.inputParameters={'cam_pixelsize_um'};
pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.name='calibrate3DsplinePSF';
pard.plugininfo.description=sprintf('Starts the gui to calibrate an experimental PSF model from bead stacks. According to: Li, Yiming, Markus Mund, Philipp Hoess, Joran Deschamps, Ulf Matti, Bianca Nijmeijer, Vilma Jimenez Sabinina, Jan Ellenberg, Ingmar Schoen, and Jonas Ries. “Real-Time 3D Single-Molecule Localization Using Experimental Point Spread Functions.” Nature Methods 15, no. 5 (April 9, 2018): 367–69. https://doi.org/10.1038/nmeth.4661.');
end
