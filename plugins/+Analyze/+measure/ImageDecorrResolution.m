classdef ImageDecorrResolution<interfaces.DialogProcessor
    %FRCresolution calculates FRC resolution
    % "Descloux, A., K. S. Grußmayer, and A. Radenovic. "Parameter-free image 
    % resolution estimation based on decorrelation analysis."
    % Nature methods (2019): 1-7."
 
    methods
        function obj=ImageDecorrResolution(varargin)           
            obj@interfaces.DialogProcessor(varargin{:}) 
            obj.inputParameters=anyRender;
            obj.showresults=true;
            obj.history=false;
        end
        
        function out=run(obj,p)
            %binning: works only with layer1.
            out=[];
            if p.takeimage
                 p.pixrec_frc=p.sr_pixrec;
            end
            locs=obj.locData.getloc({'xnm','ynm'},'layer',1,'Position','roi');
            xedges=min(locs.xnm):p.pixrec_frc:max(locs.xnm);
            yedges=min(locs.ynm):p.pixrec_frc:max(locs.ynm);
            image1=histcounts2(locs.xnm,locs.ynm,xedges,yedges);
         
            
            pps = p.pixrec_frc; % projected pixel size of 15nm
            % typical parameters for resolution estimate
            Nr = 50;
            Ng = 10;
            r = linspace(0,1,Nr);
            GPU = 1;
            image1 = apodImRect(image1,20);
            f=figure;
            figID=f.Number;
            [kcMax,A0] = getDcorr(image1,r,Ng,figID);
            ax=obj.initaxis('Decorr Resolution');
            parent=ax.Parent;
            delete(ax)
            f.Children.Parent=parent;
            delete(f)

        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end



function rs=radialsum(img)
s=size(img);
center=floor((s+1)/2);
rs=zeros(ceil(s(1)/2)+1,1);
for k=1:s(1)
    for l=1:s(2)
        d=sqrt((k-center(1)).^2+(l-center(2)).^2);
        ind=round(d)+1;
        if ind<=length(rs)
        rs(ind)=rs(ind)+img(k,l);
        end
    end
end
end
function imageo=filterimage(imagei)
sz=size(imagei);
nfac = 8;    
% Image width / Width of edge region
x=((1:sz(1))-sz(1)/2)/sz(1);
x_im = meshgrid(x);
mask = 0.5-0.5*cos(pi*nfac*x_im);          
mask(abs(x_im)<((nfac-2)/(nfac*2))) = 1;
maskt=mask.*mask';
imageo=maskt.*imagei;

end


function pard=guidef
pard.t0.object=struct('String','Image Decorr resolution, Layer 1. Uses ROI or FoV (if rendered image). ','Style','text');
pard.t0.position=[1,1];
pard.t0.Width=4;

pard.takeimage.object=struct('String','use settings from rendered image','Style','checkbox','Value',0);
pard.takeimage.position=[2,1];
pard.takeimage.Width=2;
pard.t1.object=struct('String','otherwise:','Style','text');
pard.t1.position=[3,1];
pard.t2.object=struct('String','pixelsize (nm)','Style','text');
pard.t2.position=[4,1];

pard.pixrec_frc.object=struct('String','3','Style','edit');
pard.pixrec_frc.position=[4,2];
% 

pard.plugininfo.description='FRCresolution calculates FRC resolution "Descloux, A., K. S. Grußmayer, and A. Radenovic. "Parameter-free image  resolution estimation based on decorrelation analysis.", Nature methods (2019): 1-7."';
pard.plugininfo.type='ProcessorPlugin';
end