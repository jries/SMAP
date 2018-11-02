classdef segmentNPC<interfaces.DialogProcessor&interfaces.SEProcessor
    methods
        function obj=segmentNPC(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={};
            obj.showresults=true;
        end
        
        function out=run(obj,p)  
          segmentNPCi(obj.SE,p,obj)
          out=[];
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end




function pard=guidef
pard.t1.object=struct('String','cutoff','Style','text');
pard.t1.position=[1,1];
pard.cutoff.object=struct('String','.04','Style','edit');
pard.cutoff.position=[1,2];

% pard.t2.object=struct('String','sigmaNMS','Style','text');
% pard.t2.position=[2,1];
% pard.sigmaNMS.object=struct('String','5','Style','edit');
% pard.sigmaNMS.position=[2,2];
pard.t3.object=struct('String','diameterNPC','Style','text');
pard.t3.position=[3,1];
pard.diameterNPC.object=struct('String','110','Style','edit');
pard.diameterNPC.position=[3,2];
pard.t4.object=struct('String','rim','Style','text');
pard.t4.position=[4,1];
pard.rim.object=struct('String','20','Style','edit');
pard.rim.position=[4,2];

pard.saveon.object=struct('String','saveon','Style','checkbox');
pard.saveon.position=[1,3];

pard.getmask.object=struct('String','getmask','Style','checkbox');
pard.getmask.position=[2,3];
pard.plugininfo.type='ROI_Analyze';
end

function segmentNPCi(se,p,obj)
% global se


saveon=p.saveon;
getmask=p.getmask;
%segment with free roi
% if saveon==1
if getmask
figure(24)
file=se.files(se.indexFromID(se.files,se.currentfile.ID));
imagesc(file.image.image)
colormap  hot
h = imfreehand;
bw=createMask(h);
bw(1,:)=false;bw(:,1)=false;
pf=file.image.parameters.sr_pixrec;
sf=size(file.image.image);
end
% end

%%

cutoff=p.cutoff;
sigmaNMS=3;%p.sigmaNMS;
diameterNPC=p.diameterNPC; %nm
rim=p.rim; %nm

sigmaNMS=round((sigmaNMS-1)/2)*2+1;

cells=se.cells;
pixrec=cells(1).image.parameters.sr_pixrec;

%make kernel for ring filter
rRingO=double((diameterNPC/2+rim)/pixrec);
hfilterO=fspecial('disk',rRingO);
 rRingI=double(max(1,(diameterNPC/2-rim)/pixrec));
hfilterI=fspecial('disk',rRingI);

si=size(hfilterI);
n=(si(1)-1)/2;
so=size(hfilterO);
c=(so(1)+1)/2;
hfilterIb=0*hfilterO;
hfilterIb(c-n:c+n,c-n:c+n)=hfilterI;
hfilter=hfilterO/max(hfilterO(:))-hfilterIb/max(hfilterIb(:));
hfilter=hfilter/sum(hfilter(:));

%gaussian ring kernel

n=-1*diameterNPC:pixrec:1*diameterNPC;
[X,Y]=meshgrid(n,n);

h=exp(-abs(X.^2+Y.^2-(diameterNPC/2)^2)/4/rim^2);
h=h/sum(h(:));
obj.initaxis('filter');imagesc(n,n,h);

for cn=1:length(cells)
    cell=cells(cn);
%     srim=sum(double(cell.image.image),3);
srim=0;
%     srim=cell.image.layers(1).images.rawimage.image;
    for k=1:length(cell.image.layers)
        if ~isempty(cell.image.layers(k).images)
            srim=cell.image.layers(k).images.rawimage.image+srim;
        end
    end
%     srims=size(srim);


    srimbw=sum(srim,3);
   

    % imfO=imfilter(imsrrem,hfilterO);
    % imfI=imfilter(imsrrem,hfilterI);
    imfD=imfilter(sqrt(srimbw),h);

    hfilterGauss=fspecial('gauss',21,max(1,0.7*rRingO));
    hfilterGauss2=fspecial('gauss',51,max(2,2*rRingO));
    srimf=imfilter(imfD,hfilterGauss)-imfilter(imfD,hfilterGauss2);
% srimf=imfD;
%     h=fspecial('gaussian',2*sigma,sigma);
%     srimf=filter2(h,srim);

    maximaout=NMS2DBlockCcall(srimf,sigmaNMS);
    pc=cell.image.parameters.sr_pixrec;
    sc=size(cell.image.image)/2;
    
    
    
    ysite=maximaout(:,1)*pc+cell.image.rangey(1)*1000;
    xsite=maximaout(:,2)*pc+cell.image.rangex(1)*1000;
    
    %initial mask
    if getmask
    ym=ysite-file.image.rangey(1)*1000;
    xm=xsite-file.image.rangex(1)*1000;
    xr=round(xm/pf);
    yr=round(ym/pf);
    xr(xr<1)=1;yr(yr<1)=1;xr(xr>sf(2))=1;yr(yr>sf(1))=1;
    
    linind=sub2ind(sf(1:2),yr,xr);
    indinroi=bw(linind);
    % cutoff=max(maximaout(:,3));
    indgood=maximaout(:,3)>=cutoff;
    indgood=indgood&indinroi;
    else
        indgood=maximaout(:,3)>=cutoff;
    end
    
%     newsites=maximaout(indgood,:);

xsite=xsite(indgood);
ysite=ysite(indgood);

    obj.initaxis('images')
    subplot(1,2,1)
    srimp=srim;
    maxi=myquantilefast(srimp,.99);
    srimp(srimp>maxi)=maxi;
    imagesc(srimp)
    hold on
    plot(maximaout(indgood,2),maximaout(indgood,1),'wo')
    hold off
        subplot(1,2,2)
    imagesc(srimf)
    hold on
    plot(maximaout(indgood,2),maximaout(indgood,1),'wo')
    hold off
    colorbar
    
    if saveon
    %convert to micrometers
%     poscell=cell.pos;
%     newsitesum=(newsites-srims(1)/2)*cell.globpar.pixrec/1000;
%     possites=newsitesum(:,[2 1]);
%     possites(:,1)=possites(:,1)+poscell(1);
%     possites(:,2)=possites(:,2)+poscell(2);

    for k=1:length(xsite)
        thissite=interfaces.SEsites;
        thissite.pos=[xsite(k) ysite(k) 0];
        thissite.info.cell=cell.ID;
        thissite.info.filenumber=cell.info.filenumber;
        % thissite.cellnumber=sitepar.currentcell.number;
%         thissite.number=sitepar.sitelist.cellnumber+1;
        se.addSite(thissite);

    end
    else
        waitforbuttonpress
    end
end
se.processors.preview.updateSitelist
end
