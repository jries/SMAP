classdef GaussPSF<interfaces.PSFmodel
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        locfields={'x','y','z','N'}
    end
    
    methods
        function img=render(obj,locs,xrange,yrange,pixelsizex,pixelsizey)
%             locsh.x=locs.xnm;
            p=obj.guipar;
            if length(p.analytical)==2 
                p.analytical(3)=0;
            end
            locs.sx=p.analytical(1).*sqrt(1+((locs.z+p.analytical(3))/p.analytical(2)).^2);
            locs.sy=p.analytical(1).*sqrt(1+((locs.z-p.analytical(3))/p.analytical(2)).^2);
            
            img=gaussrender_ellipt(locs,xrange+0*pixelsizex/2,yrange+0*pixelsizey/2,pixelsizex,pixelsizey);
            
        end
        function img=PSF(obj,locs)
            img=[];
%             x=locs.x;y=locs.y;z=locs.z;
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
%                 if isfield(locs,'N')
%                     output=output*locs.N;
%                 end
%                 img=output;
        end
        
        function pard=guidef(obj)
            pard.useanalytical.object=struct('String','Use analytical model [s0, zR (dz)]','Style','edit');
            pard.useanalytical.position=[1,1];
            pard.useanalytical.Width=3;
            pard.analytical.object=struct('String','100 300 300','Style','edit');
            pard.analytical.position=[1,1];
            pard.analytical.Width=3;            
            
            
            pard.calfile.object=struct('String','cal3D.mat','Style','edit');
            pard.calfile.position=[3,1];
            pard.calfile.Width=3;
            pard.load.object=struct('String','load','Style','pushbutton','Callback',{{@load_callback,obj}});
            pard.load.position=[3,4];
        end
    end
    
end

function load_callback(a,b,obj)
s=obj.getSingleGuiParameter('calfile');
[f,p]=uigetfile(s);
if f
obj.setGuiParameters(struct('calfile',[p f]));
l=load([p f]);
% obj.modelpar=(l.coeff);
obj.modelpar=permute(l.coeff,[4 1 2 3]);
end
end