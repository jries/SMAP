classdef ImposeSymmetry<interfaces.DialogProcessor&interfaces.SEProcessor
    methods
        function obj=ImposeSymmetry(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.showresults=false;
            
        end
        
        function out=run(obj,p)
            out=[];
            obj.setPar('undoModule','ImposeSymmetry');
            notify(obj.P,'backup4undo');
            
            ll=obj.locData.loc;
            if p.centerauto
                dx=median(ll.xnm);dy=median(ll.ynm);
            else
                dx=p.centerx;dy=p.centery;
            end
            if ~isfield(ll,'class')
                ll.class=0*ll.nm;
            end
            mc=max(ll.class)+1;

            x=ll.xnm-dx;
            y=ll.ynm-dy;
%             f=ll.frame;

            for k=1:p.symmetry
                a=k*2*pi/p.symmetry;
                xn=x*cos(a)-y*sin(a);
                yn=y*cos(a)+x*sin(a);
                lln=ll;
                lln.xnm=xn+dx;
                lln.ynm=yn+dy;
                lln.class=lln.class+mc*k;
                obj.locData.addLocData(lln);
            end
            obj.locData.regroup;
            obj.locData.filter;
            
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end




function pard=guidef

pard.symmetryt.object=struct('String','rotational symmetry (fold)','Style','text');
pard.symmetryt.position=[1,1];
pard.symmetryt.Width=1.5;

pard.symmetry.object=struct('String','8','Style','edit');
pard.symmetry.position=[1,2.5];
pard.symmetry.Width=.5;

pard.centert.object=struct('String','rotation center x, y (nm)','Style','text');
pard.centert.position=[2,1];
pard.centert.Width=1.5;

pard.centerx.object=struct('String','500','Style','edit');
pard.centerx.position=[2,2.5];
pard.centerx.Width=.5;
pard.centery.object=struct('String','500','Style','edit');
pard.centery.position=[2,3];
pard.centery.Width=.5;

pard.centerauto.object=struct('String','auto: median','Style','checkbox');
pard.centerauto.position=[2,3.5];
pard.centerauto.Width=1.5;

pard.plugininfo.type='ProcessorPlugin';


end