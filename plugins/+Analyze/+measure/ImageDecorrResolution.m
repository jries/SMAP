classdef ImageDecorrResolution<interfaces.DialogProcessor
% Calculates decorrelation based resolution. Careful! Returns mainly the
% size of the Gaussian kernel. Not suitable to assess image resolution in
% SMLM. According to "Descloux, A., K. S. Grußmayer, and A. Radenovic.
% "Parameter-free image  resolution estimation based on decorrelation
% analysis.", Nature methods (2019): 1-7."';

 
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
            locs=obj.locData.getloc({'xnm','ynm'},'layer',1,'Position','roi');
            xedges=min(locs.xnm):p.pixrec_frc:max(locs.xnm);
            yedges=min(locs.ynm):p.pixrec_frc:max(locs.ynm);
            switch p.reconmode.Value
                case 1 %histogram
                    image1=histcounts2(locs.xnm,locs.ynm,xedges,yedges);
                case 2 %Gaussian
                    image1=histcounts2(locs.xnm,locs.ynm,xedges,yedges);
                    h=fspecial('gaussian',ceil(7*p.gausssigma),p.gausssigma);
                    image1=imfilter(image1,h);
                case 3 %bilinear
                    [image1,pps] = smlmHist([locs.xnm-xedges(1),locs.ynm-yedges(1)],p.pixrec_frc,max(xedges(end)-xedges(1),yedges(end)-yedges(1)));
                    p.pixrec_frc=pps;
                    figure(876);imagesc(image1);
                case 4 %reconstructed image
                    image1=obj.getPar('sr_image');
                    image1=sum(image1.composite,3);
                    p.pixrec_frc=p.sr_pixrec;
            end    
            
            imclipboard('copy',(image1/max(image1(:))));

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
            
            resnm=p.pixrec_frc*2/kcMax;
            text(parent.Children,0.5,0.9,['res: ' num2str(resnm,'%2.1f') 'nm'],'FontSize',20);

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
pard.t0.object=struct('String','Image Decorr resolution, Layer 1. Uses rendered image or ROI / FoV. ','Style','text');
pard.t0.position=[1,1];
pard.t0.Width=4;

pard.t0b.object=struct('String','Careful! Returns mainly the size of the Gaussian kernel. ','Style','text');
pard.t0b.position=[2,1];
pard.t0b.Width=4;

pard.t0c.object=struct('String','Not suitable to assess image resolution in SMLM. ','Style','text');
pard.t0c.position=[3,1];
pard.t0c.Width=4;

% pard.takeimage.object=struct('String','use rendered image','Style','checkbox','Value',0);
% pard.takeimage.position=[4,1];
% pard.takeimage.Width=2;
% pard.t1.object=struct('String','otherwise:','Style','text');
% pard.t1.position=[4,3];
pard.t2.object=struct('String','pixelsize (nm)','Style','text');
pard.t2.position=[5,1];

pard.pixrec_frc.object=struct('String','3','Style','edit');
pard.pixrec_frc.position=[5,2];
% +
pard.reconmode.object=struct('String',{{'histogram','Gaussian','bilinear histogram','reconstructed image'}},'Style','popupmenu','Value',1);
pard.reconmode.position=[4,1];
pard.reconmode.Width=2;
pard.gausssigmat.object=struct('String','Sigma for Gauss','Style','text');
pard.gausssigmat.position=[6,1];
pard.gausssigmat.Width=2;
pard.gausssigma.object=struct('String','1','Style','edit');
pard.gausssigma.position=[6,3];

pard.plugininfo.description='Calculates decorrelation based resolution. Careful! Returns mainly the size of the Gaussian kernel. Not suitable to assess image resolution in SMLM. According to "Descloux, A., K. S. Grußmayer, and A. Radenovic. "Parameter-free image  resolution estimation based on decorrelation analysis.", Nature methods (2019): 1-7."';
pard.plugininfo.type='ProcessorPlugin';
end