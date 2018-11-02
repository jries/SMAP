classdef calibrateDHPSF<interfaces.DialogProcessor
%  Copyright (c)2017 Ries Lab, European Molecular Biology Laboratory,
%  Heidelberg. This file is part of Single Molecule Analysis Platform (SMAP).

%3D calibration using a double helical PSF. This plugin is based on:
%Lew et al., Easy-DHPSF open-source software for three-dimensional
%localization of single molecules with precision beyond the optical
%diffraction limit., Protocol Exchange (2013).
    
    properties
        zold
    end
    methods
        function obj=calibrateDHPSF(varargin)   
            obj@interfaces.DialogProcessor(varargin{:}) ;
             obj.inputParameters={'cam_pixelsize_nm'};
             obj.showresults=true;

        end
        function out=run(obj,p)
                    p.sigmaBounds = [1.0 1.5]*160/p.nmPerPixel;
                    p.lobeDistBounds = [3.5 10]*160/p.nmPerPixel;
                    p.boxRadius = round(7*160/p.nmPerPixel);
                    p.gaussianFilterSigma = 1.5*160/p.nmPerPixel;
                    p.minDistBetweenSMs = 7.5*160/p.nmPerPixel;
                    [outputFilePrefix, numBeads]=f_calDHPSF(obj,p);
                    obj.setPar('dhpsf_calfile',[outputFilePrefix filesep 'calibration.mat']);
                    out=outputFilePrefix;

           
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
il=imageloaderAll([p f]);
ps.nothing='';
if il.metadata.emgain>1 %metadata found, unlikely emgain = 300
    ps.offsetADU=il.metadata.offset;
    ps.EMgain=il.metadata.emgain;
    
end
if il.metadata.conversion>1
    ps.nmPerPixel=il.metadata.cam_pixelsize_um(1);
    ps.conversion=il.metadata.conversion;
end
obj.setGuiParameters(ps);
end



function pard=guidef(obj)
% pard.d3_color.object=struct('Style','checkbox','String','render in color');
% pard.d3_color.position=[2,1];
% pard.mode3D.object=struct('String',{{'astigmatic PSFx/y for MLE','astigmatic gradient fit','bi-plane'}},'Style','popupmenu');
% pard.mode3D.position=[2,3];
% pard.mode3D.Width=3;

pard.text2.object=struct('String','calibrate 3D DHPSF (easy-DHPSF)','Style','text');
pard.text2.position=[1,1];
pard.text2.Width=3;
pard.text3.object=struct('String','Conversion','Style','text');
pard.text3.position=[2,1];
pard.text4.object=struct('String','Pixel size (nm)','Style','text');
pard.text4.position=[4,1];
pard.texta.object=struct('String','EM gain','Style','text');
pard.texta.position=[3,1];
pard.textb.object=struct('String','offset (ADU)','Style','text');
pard.textb.position=[5,1];
% 

pard.conversion.object=struct('Style','edit','String','6'); 
pard.conversion.position=[2,2];
pard.conversion.Width=.5;
pard.nmPerPixel.object=struct('Style','edit','String','100'); 
pard.nmPerPixel.position=[4,2];
pard.nmPerPixel.Width=.5;
pard.EMgain.object=struct('Style','edit','String','300'); 
pard.EMgain.position=[3,2];
pard.EMgain.Width=.5;
pard.offsetADU.object=struct('Style','edit','String','100'); 
pard.offsetADU.position=[5,2];
pard.offsetADU.Width=.5;

pard.text5.object=struct('String','dz (nm)','Style','text');
pard.text5.position=[2,3];
pard.dz.object=struct('Style','edit','String','50'); 
pard.dz.position=[2,4];
pard.dz.Width=.5;

pard.text5b.object=struct('String','frames per z pos','Style','text');
pard.text5b.position=[3,3];
pard.dzrep.object=struct('Style','edit','String','1'); 
pard.dzrep.position=[3,4];
pard.dzrep.Width=.5;


pard.text6.object=struct('String','frame of 0 position','Style','text');
pard.text6.position=[4,3];


pard.framez0.object=struct('Style','edit','String','21'); 
pard.framez0.position=[4,4];
pard.framez0.Width=.5;

pard.dhpsf_beadfile.object=struct('Style','edit','String','*.tif'); 
pard.dhpsf_beadfile.position=[6,1];
pard.dhpsf_beadfile.Width=3;
pard.loadb.object=struct('Style','pushbutton','String','load','Callback',{{@load_callback,obj,'dhpsf_beadfile'}}); 
pard.loadb.position=[6,4];
% pard.reverse.object=struct('String','reverse z axis','Style','checkbox','Value',0);
% pard.reverse.position=[4,3.5];

pard.syncParameters={{'dhpsf_conversion','conversion',{'String'}},{'dhpsf_EMgain','EMgain',{'String'}},{'dhpsf_nmPerPixel','nmPerPixel',{'String'}},{'dhpsf_offsetADU','offsetADU',{'String'}}};

pard.plugininfo.name='calibrate 3D';
pard.plugininfo.description=sprintf('3D calibration using a double helical PSF. This plugin is based on: \n Lew et al., Easy-DHPSF open-source software for three-dimensional \n localization of single molecules with precision beyond the optical \n diffraction limit., Protocol Exchange (2013).');
pard.plugininfo.type='ProcessorPlugin';
end
