classdef pushLocs<interfaces.WorkflowModule
    %Inserts single-molecule localizations into workflow pipeline
    properties
        filestruc;
        locs
        preview
    end
    methods
        function obj=pushLocs(varargin)
            obj@interfaces.WorkflowModule(varargin{:});
        end
        function pard=guidef(obj)
            pard=guidef;
        end
%         function initGui(obj)
%             initGui@interfaces.WorkflowModule(obj);
%         end
        function prerun(obj,data,p)
            global intLoc2pos_ind2 intLoc2pos_locframes;
            intLoc2pos_ind2=1;
            
            p=obj.getAllParameters;
            obj.preview=obj.getPar('loc_preview');
            ind=[0 strfind(p.locfields,',') length(p.locfields)+1];
            for k=1:length(ind)-1
                locfields{k}=p.locfields(ind(k)+1:ind(k+1)-1);
            end
            obj.locData.sort('frame')
            obj.locs=obj.locData.getloc(locfields);
            obj.filestruc=obj.locData.files.file(1);
            pix_cam=obj.filestruc.info.cam_pixelsize_um*1000;
            x=double(obj.locs.xnm);
            y=double(obj.locs.ynm);
            
            
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
            elseif ~isempty(obj.locs.PSFxnm)
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
                while intLoc2pos_ind2<lf&&intLoc2pos_locframes(intLoc2pos_ind2)==frame;
                    intLoc2pos_ind2=intLoc2pos_ind2+1;
                end
                intLoc2pos_ind2=intLoc2pos_ind2-1;
                % XXXXX replace by: nm2pixLoc???
                [locpix,locpixr]=nm2pixLoc(obj.locs.xnm(ind1:intLoc2pos_ind2),obj.locs.ynm(ind1:intLoc2pos_ind2),obj.filestruc.info.cam_pixelsize_um*1000,obj.filestruc.info.roi);
                fn=fieldnames(obj.locs);
                for k=1:length(fn)
                    locout.(fn{k})=obj.locs.(fn{k})(ind1:intLoc2pos_ind2);
                end
                locout.xpix=locpix.x;
                locout.ypix=locpix.y;
               
               datout=data;%.copy;
               datout.data=locout;%.set(maxout);
               if obj.preview &&~isempty(data.data)
                     obj.setPar('preview_peakfind',locout);
               end
%                obj.output(datout); 
        end
    end
end

function [loc,locr]=nm2pixLoc(x,y,pixelsize,roi)
loc.x=(x/pixelsize(1))-roi(1);
loc.y=(y/pixelsize(end))-roi(2);
locr.x=round(loc.x);
locr.y=round(loc.y);
end


function pard=guidef
pard.locfields.object=struct('Style','edit','String','znm,xnm,ynm,frame,phot,bg');
pard.locfields.position=[1,1];
pard.locfields.Width=4;
pard.plugininfo.type='WorkflowModule'; 
pard.plugininfo.description='Inserts single-molecule localizations into workflow pipeline';
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