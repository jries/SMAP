classdef PSFmodel<interfaces.GuiModuleInterface
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        guipar
        modelpar
        X
        Y
        roisize=13;
    end
    properties(Abstract=true)
        locfields
    end
        
    
    methods (Abstract=true)
       imout=PSF(obj,locs);

    end
    methods
        function img=render(obj,locsh,xrange,yrange,varargin)
            p=parseinput(varargin{:});
            %or xrange: pixels in x, yrange pixels in y, locs already in
            %name value: pixelsizex pixelsizey lut
            %pixel units
            %z nm
            locsh=copyfields([],locsh,obj.locfields);
            if length(xrange)==1
                xrange=[0  xrange];
            end
            if length(yrange)==1
                yrange=[0  yrange];
            end            
%             if nargin>4
                sx=round((xrange(2)-xrange(1))/p.pixelsizex);
                sy=round((yrange(2)-yrange(1))/p.pixelsizey);
                xh=(locsh.x-xrange(1))/p.pixelsizex;yh=(locsh.y-yrange(1))/p.pixelsizey;
%             else
%                 sx=xrange;sy=yrange;
%                 xh=locsh.x-xrange(1); yh=locsh.y-xrange(2);
%             end
%             zh=locsh.z;
            locs1=copystructReduce(locsh,1);locs1.x=0;locs1.y=0;
            imh=obj.PSF(locs1);
            roipix=size(imh,1);
            roipixh=floor(roipix/2);

            if ~isempty(p.lut)
                col=true;
                imgh=zeros(sx+roipix,sy+roipix,3);
                if isnumeric(p.lut)
                    lut=p.lut;
                elseif ischar(p.lut)
                    lut=eval([p.lut '(256)']);
                else
                    lut=p.lut(256); %function handle
                end
                if isempty(p.zrange)
                    zrange=[-obj.modelpar.z0*obj.modelpar.dz obj.modelpar.z0*obj.modelpar.dz];
                else
                    zrange=p.zrange;
                end
                    zind=ceil((locsh.z-zrange(1))/(zrange(2)-zrange(1))*size(lut,1));
                    zind(zind<1)=1;zind(zind>size(lut,1))=size(lut,1);
                
            else
                col=false;
                imgh=zeros(sx+roipix,sy+roipix);
            end

            for k=1:length(xh)
                xr=round(xh(k)); yr=round(yh(k));
                locs1=copystructReduce(locsh,k);
                locs1.x=xh(k)-xr;locs1.y=yh(k)-yr;
                imh=obj.PSF(locs1);
                rxh=xr+1:xr+roipix;rxh=max(rxh,1);rxh=min(rxh,sx+roipix);
                ryh=yr+1:yr+roipix;ryh=max(ryh,1);ryh=min(ryh,sy+roipix);
                if col
                    for c=1:3
                        imgh(rxh,ryh,c)=imgh(rxh,ryh,c)+imh*lut(zind(k),c); 
                    end
                else
                    imgh(rxh,ryh)=imgh(rxh,ryh)+imh;    
                end
            end
            if col
                img=permute(imgh(roipixh+2:end-roipixh,roipixh+2:end-roipixh,:),[2 1 3]);
                img=img/max(img(:));
            else
                img=imgh(roipixh+2:end-roipixh,roipixh+2:end-roipixh)';
            end
        end
    end
end


function out=parseinput(varargin)

p = inputParser;
   addParameter(p,'pixelsizex',1,@isnumeric);
   addParameter(p,'pixelsizey',1,@isnumeric);
   addParameter(p,'zrange',[],@isnumeric);
   addParameter(p,'lut',[]);
%    addParameter(p,'size',[],@isnumeric);
   parse(p,varargin{:});
   out=p.Results;
end
