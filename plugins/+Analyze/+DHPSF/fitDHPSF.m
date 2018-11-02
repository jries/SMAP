classdef fitDHPSF<interfaces.DialogProcessor
%  Copyright (c)2017 Ries Lab, European Molecular Biology Laboratory,
%  Heidelberg. This file is part of Single Molecule Analysis Platform (SMAP).

%Applies 3D calibration model to fit a double helical PSF. This plugin is based on:
%Lew et al., Easy-DHPSF open-source software for three-dimensional
%localization of single molecules with precision beyond the optical
%diffraction limit., Protocol Exchange (2013).
    properties
        zold
    end
    methods
        function obj=fitDHPSF(varargin)   
            obj@interfaces.DialogProcessor(varargin{:}) ;
             obj.inputParameters={'cam_pixelsize_nm'};
             obj.showresults=true;

        end
        function out=run(obj,p)
                    
%                     [outputFilePrefix] = ...
%     f_fitSMs(dataFile,dataPath,calFile,calBeadIdx,templateFile,templateFrames,peakThreshold,...
%     darkFile,logFile,logPath,boxRadius,channel, sigmaBounds,gaussianFilterSigma,minDistBetweenSMs,...
%     lobeDistBounds,conversionGain,nmPerPixel,EMGain,templateLocs,ROI)
% p.dataFile=p.dhpsf_datafile;
% p.calFile=p.dhpsf_calfile;
% templatefile=dir([fileparts(p.calFile) filesep '*templates.mat']);
% p.templateFile=[fileparts(p.calFile) filesep templatefile(1).name];
% 
% [p.templateFrames, p.ROI, dataFile, dataPath, p.templateLocs,outputFilePrefix] = f_calSMidentification(obj,p);
 thresholdf=load(p.dhpsf_calfile);
        pfit=thresholdf.thresholdparameters;
        pfit.dataFile=p.dhpsf_datafile;
        pfit.peakThreshold=pfit.Thresholds/10000;
         pfit.peakThreshold= pfit.peakThreshold*.75;
        %remove
%         pfit.EMgain=300;
        pfit.framerange=p.framerange;
        
        
                    pfit.sigmaBounds = [1.0 1.5]*160/pfit.nmPerPixel;
                    pfit.lobeDistBounds = [3.5 10]*160/pfit.nmPerPixel;
                    pfit.boxRadius = round(7*160/pfit.nmPerPixel);
                    pfit.gaussianFilterSigma = 1.5*160/pfit.nmPerPixel;
                    pfit.minDistBetweenSMs = 7.5*160/pfit.nmPerPixel;

                    pfit.calBeadIdx=1; %later: change
        tfit=tic      ;      
        totalPSF =f_fitSMs(obj,pfit);
        timefit=toc(tfit)
        goodFit = totalPSF(:,17)>0;
%         totalPSF=totalPSF(goodFit,:);
        
        % add to loc data
        csvfile=[p.dhpsf_datafile(1:end-4) '_t' num2str(round(timefit)) '.csv']
        fid=fopen(csvfile,'w');
%                fprintf(fid, ['frame num,molecule num,x (nm),y (nm),z (nm),' ...
%             'x fid-corrected (nm),y fid-corrected (nm), z fid-corrected (nm),'...
%             'photons detected,mean background photons\n']);
                       fprintf(fid, ['frame num,molecule num,x (nm),y (nm),z (nm),' ...
            'photons detected,mean background photons,fitflag\n']);
        fclose(fid);


        
        
        dlmwrite(csvfile,totalPSF(:,[1 2 25 26 27  21 15 17]), '-append');
                    out=totalPSF;

           
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        

    end
    methods(Static)
    end
end


function load_callback(a,b,obj,name)
[f,p]=uigetfile(obj.guihandles.(name).String);
obj.guihandles.(name).String=[p f];
end

function pard=guidef(obj)
% pard.d3_color.object=struct('Style','checkbox','String','render in color');
% pard.d3_color.position=[2,1];
% pard.mode3D.object=struct('String',{{'astigmatic PSFx/y for MLE','astigmatic gradient fit','bi-plane'}},'Style','popupmenu');
% pard.mode3D.position=[2,3];
% pard.mode3D.Width=3;

pard.text2.object=struct('String','fit 3D DHPSF (easy-DHPSF)','Style','text');
pard.text2.position=[1,1];
pard.text2.Width=3;
% pard.text3.object=struct('String','Conversion','Style','text');
% pard.text3.position=[2,1];
% pard.text4.object=struct('String','Pixel size (nm)','Style','text');
% pard.text4.position=[4,1];
% pard.texta.object=struct('String','EM gain','Style','text');
% pard.texta.position=[3,1];
% pard.textb.object=struct('String','offset (ADU)','Style','text');
% pard.textb.position=[5,1];
% 

% pard.conversion.object=struct('Style','edit','String','6'); 
% pard.conversion.position=[2,2];
% pard.conversion.Width=.5;
% pard.nmPerPixel.object=struct('Style','edit','String','100'); 
% pard.nmPerPixel.position=[4,2];
% pard.nmPerPixel.Width=.5;
% pard.EMgain.object=struct('Style','edit','String','300'); 
% pard.EMgain.position=[3,2];
% pard.EMgain.Width=.5;
% pard.offsetADU.object=struct('Style','edit','String','100'); 
% pard.offsetADU.position=[5,2];
% pard.offsetADU.Width=.5;


pard.textf.object=struct('String','frames','Style','text');
pard.textf.position=[2,3];
pard.framerange.object=struct('Style','edit','String','1 Inf'); 
pard.framerange.position=[2,4];

% pard.textt.object=struct('String','Thresholds','Style','text');
% pard.textt.position=[3,3];
% pard.peakThreshold.object=struct('Style','edit','String','100'); 
% pard.peakThreshold.position=[3,4];
% pard.selectthreshold.object=struct('String','Select thresholds from data','Style','checkbox');
% pard.selectthreshold.position=[4,3];
% pard.selectthreshold.Width=2;

pard.dhpsf_datafile.object=struct('Style','edit','String','*.tif'); 
pard.dhpsf_datafile.position=[6,1];
pard.dhpsf_datafile.Width=3;
pard.loadb.object=struct('Style','pushbutton','String','load','Callback',{{@load_callback,obj,'dhpsf_datafile'}}); 
pard.loadb.position=[6,4];

pard.dhpsf_calfile.object=struct('Style','edit','String','*_thresh.mat'); 
pard.dhpsf_calfile.position=[7,1];
pard.dhpsf_calfile.Width=3;
pard.loadb2.object=struct('Style','pushbutton','String','load','Callback',{{@load_callback,obj,'dhpsf_calfile'}}); 
pard.loadb2.position=[7,4];

% pard.text5.object=struct('String','dz (nm)','Style','text');
% pard.text5.position=[2,3];
% pard.dz.object=struct('Style','edit','String','50'); 
% pard.dz.position=[2,4];
% pard.dz.Width=.5;

% pard.text5b.object=struct('String','frames per z pos','Style','text');
% pard.text5b.position=[3,3];
% pard.dzrep.object=struct('Style','edit','String','1'); 
% pard.dzrep.position=[3,4];
% pard.dzrep.Width=.5;


% pard.text6.object=struct('String','frame of 0 position','Style','text');
% pard.text6.position=[4,3];
% 
% 
% pard.framez0.object=struct('Style','edit','String','21'); 
% pard.framez0.position=[4,4];
% pard.framez0.Width=.5;

% pard.reverse.object=struct('String','reverse z axis','Style','checkbox','Value',0);
% pard.reverse.position=[4,3.5];
pard.syncParameters={{ 'dhpsf_datafile','dhpsf_datafile',{'String'}},{'dhpsf_calfile_thresh','dhpsf_calfile',{'String'}}};

pard.plugininfo.name='fit DHPSF ';
pard.plugininfo.description=sprintf('Applies 3D calibration model to fit a double helical PSF. This plugin is based on: \n Lew et al., Easy-DHPSF open-source software for three-dimensional \n localization of single molecules with precision beyond the optical \n diffraction limit., Protocol Exchange (2013).');

pard.plugininfo.type='ProcessorPlugin';
end
