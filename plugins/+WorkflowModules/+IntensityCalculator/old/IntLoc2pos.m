classdef IntLoc2pos<interfaces.WorkflowModule
    properties
        filestruc;
        locs
    end
    methods
        function obj=IntLoc2pos(varargin)
            obj@interfaces.WorkflowModule(varargin{:});
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
        end
        function prerun(obj,data,p)
            global intLoc2pos_ind2 intLoc2pos_locframes;
            intLoc2pos_ind2=1;
            
            p=obj.getAllParameters;
          
            transform=loadtransformation(obj,p.Tfile);
            obj.locs=obj.locData.getloc({'frame','xnm','ynm','znm','PSFxnm'});
            pix_cam=obj.filestruc.info.cam_pixelsize_um*1000;
            x=double(obj.locs.xnm);
            y=double(obj.locs.ynm);
            indref=transform.getRef(x,y);
            indtarget=~indref;
            obj.locs.xA=zeros(size(x));obj.locs.xB=zeros(size(x));obj.locs.yA=zeros(size(x));obj.locs.yB=zeros(size(x));
            obj.locs.xA(indref)=x(indref);obj.locs.yA(indref)=y(indref);
            [obj.locs.xA(indtarget),obj.locs.yA(indtarget)]=transform.transformCoordinatesInv(x(indtarget),y(indtarget));
            obj.locs.xB(indtarget)=x(indtarget);obj.locs.yB(indtarget)=y(indtarget);
            [obj.locs.xB(indref),obj.locs.yB(indref)]=transform.transformCoordinatesFwd(x(indref),y(indref));
            
            if ~isempty(obj.locs.znm)
                if isfield(obj.filestruc.info,'fit') %determine from fitz
                    fitz=obj.filestruc.info.fit.fitzParameters;
                    parx= [fitz(7) fitz(1) fitz(2) fitz(4) fitz(6) 0];
                    pary= [fitz(7) fitz(8) fitz(3) fitz(5) -fitz(6) 0];
                    obj.locs.PSFxpix=sigmafromz(pary,obj.locs.znm/1000,1);
                    obj.locs.PSFypix=sigmafromz(parx,obj.locs.znm/1000,1);
                else
                    d=0.42;
                    g=-0.2;
                    sx0=1.1;
                    obj.locs.PSFxpix=sigmafromz_simple(obj.locs.znm/1000,[d -g sx0]);
                    obj.locs.PSFypix=sigmafromz_simple(obj.locs.znm/1000,[d g sx0]);
                end
            else
                obj.locs.PSFxpix=double(obj.locs.PSFxnm/pix_cam(1));
                obj.locs.PSFypix=obj.locs.PSFxpix;
            end
            intLoc2pos_locframes=obj.locs.frame;
         

        end
        function datout=run(obj,data,p)
            global intLoc2pos_ind2 intLoc2pos_locframes

            
            lf=length(obj.locs.xnm);
            frame=data.frame;
                %find indices for same frame
                ind1=intLoc2pos_ind2;
                while ind1>0&&intLoc2pos_locframes(ind1)<frame && ind1<lf
                    ind1=ind1+1;
                end
                ind1=min(ind1,lf);
                intLoc2pos_ind2=ind1;
                if ~(intLoc2pos_locframes(intLoc2pos_ind2)==frame) %no localizatiaon in frame
                      datout=data;%.copy;
                     datout.data=struct('x',[]);%.set(maxout);
                    return
                end
                while intLoc2pos_ind2<=lf&&intLoc2pos_locframes(intLoc2pos_ind2)==frame
                    intLoc2pos_ind2=intLoc2pos_ind2+1;
                end
                intLoc2pos_ind2=intLoc2pos_ind2-1;
                if 0 %(intLoc2pos_ind2-ind1+1)~=(sum(obj.locs.frame==frame))
                    disp((intLoc2pos_ind2-ind1+1)-(sum(obj.locs.frame==frame)))
                    disp(frame)
                end
               [maxposA,maxposAR]=nm2pixLoc(obj.locs.xA(ind1:intLoc2pos_ind2),obj.locs.yA(ind1:intLoc2pos_ind2),obj.filestruc.info.cam_pixelsize_um*1000,obj.filestruc.info.roi);
               [maxposB,maxposBR]=nm2pixLoc(obj.locs.xB(ind1:intLoc2pos_ind2),obj.locs.yB(ind1:intLoc2pos_ind2),obj.filestruc.info.cam_pixelsize_um*1000,obj.filestruc.info.roi);  
               maxout.x=vertcat(maxposAR.x,maxposBR.x);
               maxout.y=vertcat(maxposAR.y,maxposBR.y);
               maxout.frame=frame+0*maxout.y;
               maxout.dx=vertcat(maxposA.x-maxposAR.x,maxposB.x-maxposBR.x);
               maxout.dy=vertcat(maxposA.y-maxposAR.y,maxposB.y-maxposBR.y);
               maxout.PSFxpix=vertcat(obj.locs.PSFxpix(ind1:intLoc2pos_ind2),obj.locs.PSFxpix(ind1:intLoc2pos_ind2));
               maxout.PSFypix=vertcat(obj.locs.PSFypix(ind1:intLoc2pos_ind2),obj.locs.PSFypix(ind1:intLoc2pos_ind2));
               datout=data;%.copy;
               datout.data=maxout;%.set(maxout);
%                obj.output(datout); 
                if intLoc2pos_ind2==lf
                    obj.output(datout);
                    datout.data=struct('x',[]);
                    datout.eof=true;
                end
        end
    end
end

% function [loc,locr]=nm2pixLoc(x,y,pixelsize,roi)
% loc.x=(x/pixelsize(1))-roi(1);
% loc.y=(y/pixelsize(2))-roi(2);
% locr.x=round(loc.x);
% locr.y=round(loc.y);
% end


function pard=guidef
pard.Tfile.object=struct('Style','edit','String','Tfile');
pard.Tfile.position=[1,1];
pard.Tfile.Width=1.3;
pard.pixcam.object=struct('Style','edit','String','138');
pard.pixcam.position=[2,1];
pard.PSFxnm.object=struct('Style','edit','String','138');
pard.PSFxnm.position=[2,1];
pard.plugininfo.type='WorkflowModule'; 
end

function PSFx=sigmafromz_simple(z,p)%[d g sx0]);
    PSFx=p(3).*sqrt(1+(z-p(2)).^2./p(1).^2);
end


function s=sigmafromz(par,z,B0)
par=real(par);
% parx= [d sx0 Ax Bx g mp]
s0=par(2);d=par(1);A=par(3);B=par(4)*B0;g=par(5);mp=par(6);

% s=s0*sqrt(1+(z-g+mp).^2/d^2);
s=s0*sqrt(1+(z-g+mp).^2/d^2+A*(z-g+mp).^3/d^3+B*(z-g+mp).^4/d^4);
s=real(s);
end