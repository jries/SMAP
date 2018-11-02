classdef CenterPoints<interfaces.SEEvaluationProcessor
    properties
        
    end
    methods
        function obj=CenterPoints(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function makeGui(obj)
            makeGui@interfaces.SEEvaluationProcessor(obj);
%             obj.guihandles.saveimagesb.Callback={@saveimagesb_callback,obj};
        end
        function out=run(obj,p)
            try
            out=runintern(obj,p);
            catch err
                err
                out=[];
            end
         
        end
        
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end

function pard=guidef
pard.centermode.object=struct('Style','popupmenu','String',{{'median','mean'}});
pard.centermode.position=[1,1];
pard.centermode.Width=2;

pard.centerwhat.object=struct('Style','popupmenu','String',{{'line1','line2'}});
pard.centerwhat.position=[1,3];
pard.centerwhat.Width=2;


pard.sizet.object=struct('Style','text','String','Size (nm)');
pard.sizet.position=[2,1];
pard.sizet.Width=1;

pard.size.object=struct('Style','edit','String','20');
pard.size.position=[2,2];
pard.size.Width=1;

pard.iterationst.object=struct('Style','text','String','iterations: ');
pard.iterationst.position=[3,1];
pard.iterationst.Width=1;

pard.iterations.object=struct('Style','edit','String','3');
pard.iterations.position=[3,2];
pard.iterations.Width=1;


pard.plugininfo.type='ROI_Evaluate';
pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi','layer1_','layer2_','se_sitepixelsize'};
end


function out=runintern(obj,p)

%     obj.site.pos(1)
linepos=obj.site.annotation.(p.centerwhat.selection).pos*1000;
linepos(isnan(linepos(:,1)),1)=obj.site.pos(1);
linepos(isnan(linepos(:,2)),2)=obj.site.pos(2);
lineposcorr=linepos;
locs1=obj.getLocs({'xnm','ynm','locprecznm','locprecnm'},'layer',1,'size',p.se_siteroi);
locs2=obj.getLocs({'xnm','ynm','locprecznm','locprecnm'},'layer',2,'size',p.se_siteroi);

for k=1:p.iterations
    dx1=(locs1.xnm-linepos(1,1));
    dy1=(locs1.ynm-linepos(1,2));
    dx2=(locs2.xnm-linepos(2,1));
    dy2=(locs2.ynm-linepos(2,2));
    inrange1=(dx1).^2+(dy1).^2<p.size^2;
    inrange2=(dx2).^2+(dy2).^2<p.size^2;

    switch p.centermode.selection
        case 'median'
            centerfun=@median;
        case 'mean'
            centerfun=@mean;      
    end 
    linepos(1,1)=linepos(1,1)+centerfun(dx1(inrange1));
    linepos(1,2)=linepos(1,2)+centerfun(dy1(inrange1));
    linepos(2,1)=linepos(2,1)+centerfun(dx2(inrange2));
    linepos(2,2)=linepos(2,2)+centerfun(dy2(inrange2));
    
end

switch p.centerwhat.selection
    case 'line1'
        linesel=1;
    case 'line2'
        linesel=2;
end
if ~any(isnan(linepos))
obj.site.setlinepos(linesel,linepos)
else
 obj.site.setlinepos(linesel,lineposcorr)
end
out=[];
end

