classdef makeCellGrid<interfaces.DialogProcessor&interfaces.SEProcessor
%     Generates a grid of cells for all files
    methods
        function obj=makeCellGrid(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'se_cellfov','se_sitefov'};
        end
        
        function out=run(obj,p)  
%           p=obj.getAllParameters;
            files=obj.SE.files;
            xm=min(obj.locData.loc.xnm);
            ym=min(obj.locData.loc.ynm);
            xx=max(obj.locData.loc.xnm);
            yx=max(obj.locData.loc.ynm);
            posx=xm-p.se_sitefov:p.se_cellfov:xx+p.se_sitefov;
            posy=ym-p.se_sitefov:p.se_cellfov:yx+p.se_sitefov;
            for k=1:length(files)
                for x=1:length(posx)
                    for y=1:length(posy)
                        pos=[posx(x) posy(y)]+p.se_cellfov/2;

                         currentcell=interfaces.SEsites;
                        currentcell.pos=pos;
                        currentcell.ID=0;

                        currentcell.info.filenumber=k;
                        obj.SE.addCell(currentcell);
                    end
                end
            end
            obj.SE.processors.preview.updateCelllist
                      out=[];
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end



function pard=guidef(obj)




% pard.makecells.object=struct('String','Make cell grid','Style','pushbutton','Callback',{{@makecells, obj}});
% pard.makecells.position=[1,1];
% 
% 
% % pard.segmentbeads.object=struct('String','Make cell grid','Style','pushbutton','Callback',{{@makecells, obj}});
% % pard.segmentbeads.position=[2,1];
% 
% pard.t2.object=struct('String','minimum frames','Style','text');
% pard.t2.position=[2,1];
% pard.minframes.object=struct('String','15','Style','edit');
% pard.minframes.position=[2,2];
% pard.t3.object=struct('String','diameterNPC','Style','text');
% pard.t3.position=[3,1];
% pard.diameterNPC.object=struct('String','110','Style','edit');
% pard.diameterNPC.position=[3,2];
% pard.t4.object=struct('String','rim','Style','text');
% pard.t4.position=[4,1];
% pard.rim.object=struct('String','20','Style','edit');
% pard.rim.position=[4,2];
% 
% pard.saveon.object=struct('String','saveon','Style','checkbox');
% pard.saveon.position=[1,3];
% 
% pard.getmask.object=struct('String','getmask','Style','checkbox');
% pard.getmask.position=[2,3];
pard.plugininfo.type='ROI_Analyze';
pard.plugininfo.description='Generates a grid of cells for all files';
end

