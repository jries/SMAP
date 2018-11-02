classdef segmentCellsCME<interfaces.DialogProcessor&interfaces.SEProcessor
    methods
        function obj=segmentCellsCME(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={};
        end
        
        function out=run(obj,p)  
          segmentCells(obj,p)
          out=[];
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end




function pard=guidef
pard=[];
pard.t1.object=struct('String','cutoff','Style','text');
pard.t1.position=[1,1];
pard.cutoff.object=struct('String','1','Style','edit');
pard.cutoff.position=[1,2];

pard.preview.object=struct('String','preview','Style','checkbox','Value',1);
pard.preview.position=[3,1];
pard.plugininfo.type='ROI_Analyze';
% 
% pard.t2.object=struct('String','sigmaNMS','Style','text');
% pard.t2.position=[2,1];
% pard.sigmaNMS.object=struct('String','5','Style','edit');
% pard.sigmaNMS.position=[2,2];
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
end

function segmentCells(obj,p)
se=obj.SE;
locData=obj.locData;
files=se.files;
hgauss=fspecial('gaussian',15,1);
hgauss2=fspecial('gaussian',45,5);
cutoff=2;
cutoffmax=p.cutoff;
    hf=initaxis(p.resultstabgroup,'find cells');
    
if p.preview
    rangef=se.currentfile.ID;
else
rangef=1:length(files);
end

 for kf=rangef
    pr=files(kf).info.cam_pixelsize_um*1000;
    locs=locData.getloc({'xnm','ynm'},'filenumber',files(kf).ID,'layer',1,'position','all');
    roi([1 3])=files(kf).info.roi([1 3])*pr(1);
    roi([2 4])=files(kf).info.roi([2 4])*pr(2);
    rx=[roi(1) roi(1)+roi(3)];
    ry=[roi(2) roi(2)+roi(4)];
    imh=(myhist2(locs.xnm,locs.ynm,pr(1),pr(2), rx,ry))';
    imfs=imfilter(sqrt(imh),hgauss);
    imf=imfilter((imh),hgauss);
%     imbw=imf>cutoff;
%     imf(imf<cutoff)=0;
    imf2s=imfilter(double(imfs),hgauss2);
    imf2=imfilter(double(imf),hgauss2);
%     imf2=imfilter(double(imbw),hgauss2);
    
    maxima=maximumfindcall(imf2s);
    for k=1:length(maxima)
        maxima(k,3)=imf2(maxima(k,2),maxima(k,1));
    end
    indco=maxima(:,3)>cutoffmax;
    
    cx=maxima(indco,1)*pr(1)+rx(1);
    cy=maxima(indco,2)*pr(2)+ry(1);
    s=size(imh);
    implot=zeros(s(1),s(2),3);
    implot(:,:,1)=imf2/max(imf2(:));
    implot(:,:,2)=imf/max(imf(:))*3;
    imagesc(rx,ry,implot,'Parent',hf)
    hold on
    plot(cx,cy,'m+','Parent',hf);
    hold off
%     waitforbuttonpress
    if ~p.preview
    for kc=1:length(cx)
    
        currentcell=interfaces.SEsites;
        currentcell.pos=[cx(kc),cy(kc)];
        currentcell.ID=0;

        currentcell.info.filenumber=kf;
        obj.SE.addCell(currentcell);
    end
    end
    
    
end

se.processors.preview.updateCelllist
end
