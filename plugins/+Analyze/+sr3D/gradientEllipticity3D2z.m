classdef gradientEllipticity3D2z<interfaces.DialogProcessor
    properties
        cal3D
        outsx2sy2
    end
    methods
        function obj=gradientEllipticity3D2z(varargin)        
           obj@interfaces.DialogProcessor(varargin{:}) ;
             obj.inputParameters={'cam_pixelsize_nm'};
        end
        function initGui(obj)
            set(obj.guihandles.savebutton,'Callback',{@loadcalfile,obj})
        end
        function out=run(obj,p)
            out=[];
            locs=obj.locData.getloc({'gradient3Dellipticity'});
            if isempty(locs.gradient3Dellipticity)
                error('gradient3Dellipticity required')
            end
            epsl=log(locs.gradient3Dellipticity);
            
            switch p.fitmode.selection
                case 'Linear fit'
                    b=p.fitpol;
                    z=b(1)+b(2)*epsl;

                otherwise
                    l=load(p.calfile);
                    fitp=l.outforfit.polynomial;
                    z=fitp(epsl);
                    
                    
            end
            obj.locData.loc.znm=z*1000*p.refractiveIndexMismatch;
            obj.locData.regroup;
            obj.setPar('locFields',fieldnames(obj.locData.loc))
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end


function loadcalfile(a,b,obj)
fn=obj.guihandles.calfile.String;
[f,p]=uigetfile(fn);
if f
l=load([p f]);
txt=sprintf([num2str(l.outforfit.linear(1)) '\t' num2str(l.outforfit.linear(2))]); 
obj.setPar('fit_gradient3Dellipticity',(txt))
obj.guihandles.calfile.String=[p f];
end
end

function pard=guidef

pard.fitmode.object=struct('String',{{'Polynomial fit from file', 'Linear fit'}},'Style','popupmenu');
pard.fitmode.position=[1,1];
pard.fitmode.Width=2;

pard.pol_check.object=struct('String','linear fit:','Style','text');
pard.pol_check.position=[2,1];
pard.pol_check.Width=1;


pard.fitpol.object=struct('Style','edit','String','-.055 .35'); 
pard.fitpol.position=[2,2];
pard.fitpol.Width=2;

pard.t1.object=struct('Style','text','String','refractive Index Mismatch factor (<=1)'); 
pard.t1.position=[3,1];
pard.t1.Width=2;
pard.refractiveIndexMismatch.object=struct('Style','edit','String','.8'); 
pard.refractiveIndexMismatch.position=[3,3];

pard.calfile.object=struct('Style','edit','String','settings/cal_3DAcal.mat');
pard.calfile.position=[4,1];
pard.calfile.Width=3;

pard.savebutton.object=struct('String','load','Style','pushbutton');
pard.savebutton.position=[4,4];

pard.syncParameters={{'fit_gradient3Dellipticity','fitpol',{'String'}}};

pard.plugininfo.name='z from gradient3D fit';
pard.plugininfo.type='ProcessorPlugin';
end
