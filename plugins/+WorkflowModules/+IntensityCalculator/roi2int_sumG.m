classdef roi2int_sumG<interfaces.GuiModuleInterface 
    %determines intensity and background in a ROI around localizations by
    %summing up the ROI and a larger ROI;
    methods
        function obj=roi2int_sumG(varargin)
            obj@interfaces.GuiModuleInterface(varargin{:});
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function out=evaluate(obj,p,img,info)
            if ~isempty(info.bgim)
                bg=info.bgim;
            else 
                bg=info.bg;
            end
            out=roi2int_sum(obj,p,img,bg);
        end
        function prerun(obj,p)
        end
    end
end

function out=roi2int_sum(obj,p,roi,bg)
sim=size(roi);
mp=round(sim(1:2)-1)/2;
dn=round((p.roisize_sum+1)/2);
if length(sim)==2
    sim(3)=1;
end
out=zeros(sim(3),2,'single');

bgnorm=(2*dn+1)^2;
for k=1:sim(3)
    imh=roi(mp(1)-dn:mp(1)+dn,mp(2)-dn:mp(2)+dn,k);
    imsum=sum(imh(:));
    bgsum=bg(k)*bgnorm;
    out(k,1)=imsum-bgsum;
    out(k,2)=bg(k);
end
end

function pard=guidef(obj)
pard.t1.object=struct('Style','text','String','roisize');
pard.t1.position=[1,1];
pard.t1.TooltipString='Roi size around localizations for summing';
% pard.t1.Width=0.5;
pard.roisize_sum.object=struct('Style','edit','String','3');
pard.roisize_sum.position=[1,2];
pard.roisize_sum.TooltipString=pard.t1.TooltipString;
info.prefix='sum';
info.name='sum';
info.fields={'sum_n','sum_bg'};
pard.plugininfo=info;
pard.plugininfo.type='WorkflowIntensity'; 
end
