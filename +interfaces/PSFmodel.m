classdef PSFmodel<interfaces.GuiModuleInterface
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        guipar
        modelpar
        X
        Y
    end
    properties(Abstract=true)
        locfields
    end
        
    
    methods (Abstract=true)
       imout=PSF(obj,locs);

    end
    methods
        function img=render(obj,locsh,xrange,yrange,pixelsizex,pixelsizey)
            locsh=copyfields([],locsh,obj.locfields);
            sx=round((xrange(2)-xrange(1))/pixelsizex);
            sy=round((yrange(2)-yrange(1))/pixelsizey);

            xh=(locsh.x-xrange(1))/pixelsizex;yh=(locsh.y-yrange(1))/pixelsizey;
%             zh=locsh.z;
            locs1=copystructReduce(locsh,1);locs1.x=0;locs1.y=0;
            imh=obj.PSF(locs1);
            roipix=size(imh,1);
            roipixh=floor(roipix/2);


            imgh=zeros(sx+roipix,sy+roipix);

            for k=1:length(xh)
                xr=round(xh(k)); yr=round(yh(k));
                locs1=copystructReduce(locsh,k);
                locs1.x=xh(k)-xr;locs1.y=yh(k)-yr;
                imh=obj.PSF(locs1);
                rxh=xr+1:xr+roipix;rxh=max(rxh,1);rxh=min(rxh,sx+roipix);
                ryh=yr+1:yr+roipix;ryh=max(ryh,1);ryh=min(ryh,sy+roipix);
                imgh(rxh,ryh)=imgh(rxh,ryh)+imh;      
            end
            img=imgh(roipixh+2:end-roipixh,roipixh+2:end-roipixh)';
        end
    end
end

