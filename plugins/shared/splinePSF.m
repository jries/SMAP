classdef splinePSF<interfaces.PSFmodel
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        locfields={'x','y','z','N'}
    end
    
    methods
        function img=PSF(obj,locs)
            x=locs.x;y=locs.y;z=locs.z;
           Npixels = 13;  
           %convert to pixel unit, center = 0
                x=x+Npixels/2-1;y=y+Npixels/2-1;
                %
                dz=50;slicez=25;
                z=(z/dz+slicez/2);
                
                xc = -2*(x - Npixels/2+0.5);
                yc = -2*(y - Npixels/2+0.5);
                zc = z -floor(z);
               
                
                output = (zeros(Npixels,Npixels));
                
                coeff=obj.modelpar;
%                 coeffp=permute(coeff,[4 1 2 3]);
                spline_xsize = size(coeff,2);
                spline_ysize = size(coeff,3);
                spline_zsize = size(coeff,4);
                off = ((spline_xsize+1)-2*Npixels)/2;
                
                xstart=floor(xc);xc=xc-xstart;
                ystart=floor(yc);yc=yc-ystart;
                zstart = floor(z);
                 [delta_f]=computeDelta3Dj(single(xc),single(yc),single(zc));
                 
                 
                for ii = 0:Npixels-1
                    for jj = 0:Npixels-1
                         model = fAt3Dj(2*ii+xstart+off,2*jj+ystart+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f,coeff);                
                         output(ii+1,jj+1)=model;
                    end
                end
                
                
                if isfield(locs,'N')
                    output=output*locs.N/sum(output(:));
                end
                img=output;
                
        end
        
        function pard=guidef(obj)
            
            pard.calfile.object=struct('String','cal3D.mat','Style','edit');
            pard.calfile.position=[1,1];
            pard.calfile.Width=3;
            pard.load.object=struct('String','load','Style','pushbutton','Callback',{{@load_callback,obj}});
            pard.load.position=[1,4];
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