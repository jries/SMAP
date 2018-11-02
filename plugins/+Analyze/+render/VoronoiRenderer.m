classdef VoronoiRenderer<interfaces.Renderer
    methods
        function obj=VoronoiRenderer(varargin)        
            obj@interfaces.Renderer(varargin{:}) ;
            obj.inputParameters={'sr_pixrec','numberOfLayers','sr_pos','sr_size','layers','sr_layerson'};
        end
        
        function out=renderfunction(obj,locs,p)
             pos=[p.sr_pos(1)-p.sr_size(1) p.sr_pos(2)-p.sr_size(2) 2*p.sr_size(1) 2*p.sr_size(2)];
                out=vrender(locs,p.sr_pixrec,pos);
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end
function imageo=vrender(locs,px,pos)
                warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId')
                A=zeros(length(locs.xnm),9);
                A(:,2)=locs.frame;
                A(:,4)=locs.xnm-pos(1);
                A(:,5)=locs.ynm-pos(2);  
                A(:,7)=locs.phot;
                if ~isempty(locs.znm)
                    A(:,6)=locs.znm;
                end
                [fov,indmax]=max(pos(3:4));
                
                
                imvt=drawVor(A,px,fov);
                %reduce size;
                ds=round((pos(4)-pos(3))/px/2);
%                 mod((pos(4)-pos(3))/px,2)
                if ds>0
                    imvt(:,1:ds)=[];
                    imvt(:,end-ds+1:end)=[];
                else
                    imvt(1:ds,:)=[];
                    imvt(end-ds:end,:)=[];                    
                end
                    
                
                imageo.image=imvt;
%                 imageo.lut=lutall;
                imageo.rangex=[pos(1) pos(1)+fov];
                imageo.rangey=[pos(2) pos(2)+fov];
                imageo.numberOfLocs=length(locs.xnm);
                imageo.istiff=false;
end
function pard=guidef
pard.t1.object=struct('String','Voronoi Renderer from SharpViSu','Style','text');
pard.t1.position=[1,2];
pard.t1.Width=4;
pard.plugininfo.type='Renderer';

end

