classdef Sxsy2zcomplete<interfaces.DialogProcessor
    properties
%         cal3D
%         outsx2sy2
%         spline
    end
    methods
        function obj=Sxsy2zcomplete(varargin)        
           obj@interfaces.DialogProcessor(varargin{:}) ;
             obj.inputParameters={'cam_pixelsize_nm'};
           obj.showresults=false;
            obj.history=true;
        end
        function initGui(obj)
            set(obj.guihandles.savebutton,'Callback',{@loadcalfile,obj})
        end
        function out=run(obj,p)
            out=[];
            locs=obj.locData.getloc({'PSFxnm','PSFynm','xnm','ynm'});
            if isempty(locs.PSFynm)
                error('PSFynm required')
            end
            
            ld=load(p.calfile);
            SXY=ld.SXY;
            s=size(SXY);
            % determine Z for depth-dependent calibration from set focus
            Z=1;
            zall=zeros(size(locs.xnm),'single');
            dall=zall+inf;
            zall=zall-inf; %default:
            dzall=zall;
            for X=1:s(1)
                for Y=1:s(2)
                    SXYh=SXY(X,Y);
                    if SXYh.Xrange(1)==SXYh.Xrangeall(1), SXYh.Xrange(1)=-inf;end
                    if SXYh.Xrange(end)==SXYh.Xrangeall(end), SXYh.Xrange(end)=inf;end
                    if SXYh.Yrange(1)==SXYh.Yrangeall(1), SXYh.Yrange(1)=-inf;end
                    if SXYh.Yrange(end)==SXYh.Yrangeall(end), SXYh.Yrange(end)=inf;end
                    indh=locs.xnm>SXYh.Xrange(1)&locs.xnm<SXYh.Xrange(2)&locs.ynm>SXYh.Yrange(1)&locs.ynm<SXYh.Yrange(2);
                    sxpix=locs.PSFxnm(indh)./p.cam_pixelsize_nm;
                    sypix=locs.PSFynm(indh)./p.cam_pixelsize_nm;
                    switch p.zcalibmode.Value
                        case 1 %spline lut
                            [zh,d]=zfromSXSYLut(SXYh.splineLUT,sxpix,sypix);
                            d2=d/sqrt(2);
                            [zh2]=horzcat(zh,zfromSXSYLut(SXYh.splineLUT,sxpix+d2,sypix+d2),zfromSXSYLut(SXYh.splineLUT,sxpix-d2,sypix-d2),...
                                zfromSXSYLut(SXYh.splineLUT,sxpix+d2,sypix-d2),zfromSXSYLut(SXYh.splineLUT,sxpix-d2,sypix+d2));
                            dzall(indh)=std(zh2,[],2);
                            
                            dall(indh)=d;
                        case 2 %sx2-sy2
                            zh=zfromSx2_Ss2(SXYh.Sx2_Sy2,sxpix,sypix);
                            
                    end
                    
                    if p.reversez
                        zh=-zh;
                    end
                    zall(indh)=zh; 
                end
            end
            
            switch p.zcalibmode.Value
                case 1 %spline lut
                    obj.locData.loc.zfromPSFdistance=single(dall);
                    locprecnmz=sqrt((obj.locData.loc.locprecnm*2.5).^2+dzall.^2/4);
                case 2 %sx2-sy2
                    locprecnmz=obj.locData.loc.locprecnm*2.5;
            end
            
            znm=zall*p.refractiveIndexMismatch;
            
            if p.setzt
                znm(isinf(znm))=p.setz;
            end
            if p.setszt
                locprecnmz(isinf(locprecnmz))=p.setsz;
            end
          
            obj.locData.loc.znm=single(znm);
            obj.locData.loc.locprecnmz=single(locprecnmz);
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
    if ~isfield(l,'outforfit') && ~isfield(l,'SXY')
        msgbox('no 3D data recognized. Select other file.');
    end
%     load([p f])
%     if exist('cal3D','var')
%         obj.cal3D=cal3D;
%     end
%     if exist('outspline','var')
%         obj.spline=outspline;
%     end
%     obj.outsx2sy2=outsx2sy2;
    obj.guihandles.calfile.String=[p f];
end
end

function pard=guidef
% pard.d3_color.object=struct('Style','checkbox','String','render in color');
% pard.d3_color.position=[2,1];
pard.t1.object=struct('Style','text','String','refractive Index Mismatch factor (<=1)'); 
pard.t1.position=[1,1];
pard.t1.Width=2;
pard.refractiveIndexMismatch.object=struct('Style','edit','String','.8'); 
pard.refractiveIndexMismatch.position=[1,3];

pard.reversez.object=struct('String','reverse z-direction (M!)','Style','checkbox');
pard.reversez.position=[2,1];
pard.reversez.Width=2;

pard.zcalibmode.object=struct('String',{{'spline look up table','Sx2-Sy2'}},'Style','popupmenu');
pard.zcalibmode.position=[3,1];
pard.zcalibmode.Width=2;


pard.zposobjective_check.object=struct('String','nominal z-position (um) of focus  (gel calib)','Style','text');
pard.zposobjective_check.position=[4,1];
pard.zposobjective_check.Width=2;

pard.zposobjective.object=struct('Style','edit','String','1000'); 
pard.zposobjective.position=[4,3];
% 
% pard.pol_check.object=struct('String','give polynome here:','Style','checkbox');
% pard.pol_check.position=[2,1];
% pard.pol_check.Width=2;

% 
% pard.fitpol.object=struct('Style','edit','String','1 0 0'); 
% pard.fitpol.position=[2,3];

pard.setzt.object=struct('Style','checkbox','String','set unassigned znm to','Value',1); 
pard.setzt.position=[6,1];
pard.setzt.Width=1.5;
pard.setz.object=struct('Style','edit','String','-2000'); 
pard.setz.position=[6,2];
pard.setz.Width=.5;
pard.setszt.object=struct('Style','checkbox','String','set unassigned locprec z nm to','Value',1); 
pard.setszt.position=[6,3];
pard.setszt.Width=1.5;
pard.setsz.object=struct('Style','edit','String','300'); 
pard.setsz.position=[6,4];
pard.setsz.Width=0.5;

pard.calfile.object=struct('Style','edit','String','settings/cal_3DAcal.mat');
pard.calfile.position=[8,1];
pard.calfile.Width=3;


pard.savebutton.object=struct('String','load','Style','pushbutton');
pard.savebutton.position=[8,4];

% pard.usespline.object=struct('String','Use spline approximation','Style','checkbox');
% pard.usespline.position=[6,1];
% pard.usespline.Width=2;



pard.syncParameters={{'cal3D_file','calfile',{'String'}}};

pard.plugininfo.type='ProcessorPlugin';

end
