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
        function out=evaluate(obj,roi1,bg1)
            out=roi2int_sum(obj,roi1,bg1);
        end
        function prerun(obj,p)
            global roi2int_sumG_parameters;
            roi2int_sumG_parameters=obj.getAllParameters;
            roi2int_sumG_parameters.initialize=true;
        end
    end
end

function out=roi2int_sum(obj,roi,bg)
global roi2int_sumG_parameters;
if roi2int_sumG_parameters.initialize
%     roi2int_sumG_parameters=obj.getAllParameters;
    roisize=roi2int_sumG_parameters.roisize_sum;
    sim=size(roi);
    roi2int_sumG_parameters.mp=round(sim(1:2)-1)/2;
    roi2int_sumG_parameters.dn=min(min(roi2int_sumG_parameters.mp)+1,round((roisize+1)/2));
    roi2int_sumG_parameters.initialize=false;
end
mp=roi2int_sumG_parameters.mp;
dn=roi2int_sumG_parameters.dn;

sim=size(roi);
if length(sim)==2
    sim(3)=1;
end
out=zeros(sim(3),2,'single');

bgnorm=(2*dn+1)^2;
for k=1:sim(3)
    imh=roi(mp(1)-dn:mp(1)+dn,mp(2)-dn:mp(2)+dn,k);
    bgh=bg(mp(1)-dn:mp(1)+dn,mp(2)-dn:mp(2)+dn,k);
    imsum=sum(imh(:));
    bgsum=sum(bgh(:));
%     imsum=squeeze(sum(sum(roi(mp(1)-dn:mp(1)+dn,mp(2)-dn:mp(2)+dn,k),1),2));
%     bgsum=squeeze(sum(sum(bg(mp(1)-dn:mp(1)+dn,mp(2)-dn:mp(2)+dn,k),1),2));
    out(k,1)=imsum-bgsum;
    out(k,2)=bgsum/bgnorm;
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
