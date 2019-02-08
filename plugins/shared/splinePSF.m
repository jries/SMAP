classdef splinePSF<interfaces.PSFmodel
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        locfields={'x','y','z','N','bg'}
    end
    
    methods
        function img=PSF(obj,locs)
            roisize=min(obj.roisize,size(obj.modelpar.coeff,1))+2;
            dn=round((roisize-1)/2);
            if isstruct(locs)
                if ~isfield(locs,'N')
                     N=ones(size(locs.x));
                else
                    N=locs.N;
                end
                if ~isfield(locs,'bg')
                     bg=zeros(size(locs.x));
                else
                    bg=locs.bg;
                end
                cor=[locs.x+dn,locs.y+dn,-locs.z/obj.modelpar.dz+obj.modelpar.z0];
            else
                cor=zeros(size(locs),'single');
                cor(:,1:2)=locs(:,1:2)+dn;
                cor(:,3)=-locs(:,3)/obj.modelpar.dz+obj.modelpar.z0;
                if size(locs,2)>3
                    N=locs(:,4);
                else
                    N=1;
                end
                if size(locs,2)>4
                    bg=locs(:,5);
                else
                    bg=0;
                end
            end
            
            imgi=simSplinePSF_call(roisize,obj.modelpar.coeff,N,bg,cor);
            img=imgi(2:end-1,2:end-1,:);
            
%            Npixels = 13;  
%            %convert to pixel unit, center = 0
%                 x=x+Npixels/2-1;y=y+Npixels/2-1;
%                 %
%                 dz=50;slicez=25;
%                 z=(z/dz+slicez/2);
%                 
%                 xc = -2*(x - Npixels/2+0.5);
%                 yc = -2*(y - Npixels/2+0.5);
%                 zc = z -floor(z);
%                
%                 
%                 output = (zeros(Npixels,Npixels));
%                 
%                 coeff=obj.modelpar;
% %                 coeffp=permute(coeff,[4 1 2 3]);
%                 spline_xsize = size(coeff,2);
%                 spline_ysize = size(coeff,3);
%                 spline_zsize = size(coeff,4);
%                 off = ((spline_xsize+1)-2*Npixels)/2;
%                 
%                 xstart=floor(xc);xc=xc-xstart;
%                 ystart=floor(yc);yc=yc-ystart;
%                 zstart = floor(z);
%                  [delta_f]=computeDelta3Dj(single(xc),single(yc),single(zc));
%                  
%                  
%                 for ii = 0:Npixels-1
%                     for jj = 0:Npixels-1
%                          model = fAt3Dj(2*ii+xstart+off,2*jj+ystart+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f,coeff);                
%                          output(ii+1,jj+1)=model;
%                     end
%                 end
%                 
%                 
%                 if isfield(locs,'N')
%                     output=output*locs.N/sum(output(:));
%                 end
%                 img=output;
                
        end
        
        function pard=guidef(obj)
            
            pard.calfile.object=struct('String','cal3D.mat','Style','edit');
            pard.calfile.position=[1,1];
            pard.calfile.Width=3;
            pard.load.object=struct('String','load','Style','pushbutton','Callback',{{@load_callback,obj}});
            pard.load.position=[1,4];
        end
        function loadmodel(obj,file)
            l=load(file);
            obj.modelpar.coeff=single(l.SXY(1).cspline.coeff{1});
            obj.modelpar.dz=l.SXY(1).cspline.dz;
            obj.modelpar.z0=l.SXY(1).cspline.z0;
            obj.modelpar.x0=l.SXY(1).cspline.x0;
        end
        function crlb=crlb(obj,N,bg,z,rois)
            if nargin<5
                rois=obj.roisize;
            end
            coeff=obj.modelpar.coeff; 
            x=rois/2;y=rois/2;
            v1=ones(length(z),1);
            zh=-(z/obj.modelpar.dz)+obj.modelpar.z0;
            coords=[v1*x , v1*y , N, bg, zh];
            crlb=CalSplineCRLB(coeff, rois, coords);
            crlb(:,5)=crlb(:,5)*obj.modelpar.dz.^2;
        end
    end
    
end

function load_callback(a,b,obj)
s=obj.getSingleGuiParameter('calfile');
[f,p]=uigetfile(s);
if f
    obj.loadmodel([p f]);
    obj.setGuiParameters(struct('calfile',[p f]));
% obj.modelpar=permute(l.coeff,[4 1 2 3]);
end
end