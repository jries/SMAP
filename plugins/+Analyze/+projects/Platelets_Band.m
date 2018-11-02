classdef Platelets_Band<interfaces.DialogProcessor
    methods
        function obj=Platelets_Band(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.inputParameters={'layer1_','sr_pixrec','sr_layerson'};

        end
        
        function out=run(obj,p)
          
%             notify(obj.locData,'updateParameters');
%             obj.setParameters(obj.locData.parameters);
%             obj.setParameters(obj.locData.layer(1).parameters.rec_addpar)
%             p=copyfields(obj.getGuiParameters.par,p);
            
%             p=copyfields(pl,p);
            [locs1,~,hroi]=obj.locData.getloc({'xnm','ynm','znm','locprecznm','locprecnm'},'layer',1,'position','roi');
            p.posRoi=hroi.getPosition*1000;
            p.file=obj.locData.files.file(1).name;
            if 0%p.sr_layerson(2)
%                 locs2=obj.locData.getlocRoi('xnm','ynm','znm','locprecznm','locprecnm',2);
            else
                locs2=[];
            end
            reconstruct_sideband(locs1,locs2,p)
%             locs2=obj.locData.getlocRoi('xnm','ynm','znm','locprecznm','locprecnm',2);
            
%             sideview_reconstruct(locs,p);
             
            out=0;
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        

    end
end

function reconstruct_sideband(locs1,locs2,p)
locs1=getPol(locs1,p);
dtheta=2*pi/p.numberOfThetas;
thetas=-pi:dtheta:pi+p.numoverlap*dtheta;
p.range=[0 (p.posRoi(4)-p.posRoi(3))/2];
rangez=[p.zmin p.zmax];
ax1=initaxis(p.resultstabgroup,'sideview');

locs1.xnmline=locs1.rho;
srall=sideview_reconstruct(locs1,p);
sim=size(srall);


pos.x=locs1.xnm;
pos.y=locs1.ynm;
pos.sx=locs1.locprecnm;
pos.sy=locs1.locprecnm;

   midp(2)=(p.posRoi(3)+p.posRoi(4))/2;
   midp(1)=(p.posRoi(1)+p.posRoi(2))/2;
    
rangex=sort([p.posRoi(1) p.posRoi(2)]);
rangey=sort([p.posRoi(3) p.posRoi(4)]);


syim=gaussrender_ellipt(pos,rangex, rangey, p.sr_pixrec, p.sr_pixrec);
syimfit=imresize(syim,[sim(1),NaN]);
syimfit=syimfit/myquantile(syimfit(:),0.995)*max(srall(:))/p.numberOfThetas*p.numoverlap*2;
sxy=size(syimfit);
allimages=zeros(sim(1),sim(2)+sxy(2)+1,p.numberOfThetas);



for k=1:p.numberOfThetas
    thetas(k)
%     indg=locs1.theta>thetas(k)&locs1.theta<thetas(k+1+1);
    indg=filtertheta(locs1.theta,thetas(k),thetas(k+p.numoverlap));
    locs1f=filterLocs(locs1,indg);
    locs1f.xnmline=locs1f.rho;
    
    srimxz=sideview_reconstruct(locs1f,p);

    
    syimfitl=addthetaline(syimfit,thetas(k+round(p.numoverlap/2)));
%     figure(33)
%     imagesc(syimfitl)
%     waitforbuttonpress
    
    collage=horzcat(syimfitl,ones(size(srimxz,1),1,'single'),srimxz);
    allimages(:,:,k)=collage;
    axes(ax1)
    imagesc(collage);
    axis equal
    colormap hot
    drawnow
end

    if 0%p.d3_color
        options.color=true;
    else
        options.color=false;
%         imstack=squeeze(sum(imstack,3));
    end
    options.message=true;
    options.comp='lzw';

    [p,f]=fileparts(p.file);

    [f,p]=uiputfile([p filesep f '_platelet.tif']);
    fileout=[p f];
    if f
    imout=uint8(allimages/max(allimages(:))*(2^8-1));
    saveastiff(imout,fileout,options)
    end

end

function ind=filtertheta(thetas,theta1,theta2)
 ind=thetas>theta1&thetas<theta2;
if theta2>pi
    ind=ind|thetas<theta2-2*pi;
end
end

function imout=addthetaline(imin,theta)
    s=size(imin);
    maxim=max(imin(:));
    imout=imin;
    
    r=0:sqrt(s(1)^2+s(2)^2);
    [x,y]=pol2cart(theta,r);
    i=round(x+s(1)/2);
    j=round(y+s(2)/2);
%     if theta<pi/2
%         xall=-s(1)/2+1:0;
%         sign=1;
%     elseif theta<0
%          xall=0:s(1);
%         sign=1;
%     elseif theta<pi
%         xall=0:s(1);
%         sign=1;
%     else
%         xall=0:s(1);
%         sign=1;
%     end
%        
        
    for k=1:length(i)
        if i(k)<=s(2)&&j(k)<=s(1)&&i(k)>0&&j(k)>0
        imout(j(k),i(k))=maxim;
        end
    end
    
end

function locso=filterLocs(locs,indg)
locso=locs;
fn=fieldnames(locs);
for k=1:length(fn)
    locso.(fn{k})=locs.(fn{k})(indg);
end
end

function locs=getPol(locs,p)
    midp(2)=(p.posRoi(3)+p.posRoi(4))/2;
    midp(1)=(p.posRoi(1)+p.posRoi(2))/2;
    
    x=(locs.xnm-midp(1));
    y=(locs.ynm-midp(2));
    [theta,rho]=cart2pol(x,y);
    locs.theta=theta;
    locs.rho=rho;
    
end

function srimxz=sideview_reconstruct(locs,p)
gf=p.layer1_.gaussfac;
% figure(222)
pos.x=locs.xnmline;
pos.y=locs.znm;

rangex=p.range;
rangey=[p.zmin p.zmax];

if p.pixauto
pixelsx=p.sr_pixrec;pixelsy=p.sr_pixrec;
else
    pixelsx=p.pixrecset;pixelsy=p.pixrecset;
end
    
ming=max(p.layer1_.mingausspix*pixelsx,p.layer1_.mingaussnm);
pos.sx=max(locs.locprecnm*gf,ming);
pos.sy=max(locs.locprecznm*gf,ming);

srimxz=gaussrender_ellipt(pos,rangex, rangey, pixelsx, pixelsy);%, lut,rangec,template);
% recgui.initaxis(p.resultstabgroup,'side view')
% imagesc(rangex,rangey,double(srimxz));
% axis equal  image xy
% colormap hot
% recgui.initaxis(p.resultstabgroup,'top view')
% pos.y=locs.ynmline;
% pos.sy=max(locs.locprecnm*gf,ming);
% rangey=[min(pos.y) max(pos.y)];
% srimxy=gaussrender_ellipt(pos,rangex, rangey, pixelsx, pixelsy);%, lut,rangec,template);
% imagesc(rangex,rangey,double(srimxy));
% axis equal image
% axis xy
% colormap hot
% plot(locs.xnm,locs.ynm,'.')
end


function pard=guidef
pard.text1.object=struct('String','parameters','Style','text');
pard.text1.position=[1,1];

pard.text2.object=struct('String','zmin','Style','text');
pard.text2.position=[2,1];
pard.text3.object=struct('String','zmax','Style','text');
pard.text3.position=[3,1];

pard.zmin.object=struct('Style','edit','String',-400); 
pard.zmin.position=[2,2];
pard.zmax.object=struct('Style','edit','String',400); 
pard.zmax.position=[3,2];

pard.pixauto.object=struct('Style','checkbox','String','pixelsize auto','Value',1);
pard.pixauto.position=[4,1];
pard.pixrecset.object=struct('Style','edit','String',5); 
pard.pixrecset.position=[4,2];

pard.text4.object=struct('String','number of frames','Style','text');
pard.text4.position=[5,1];
pard.numberOfThetas.object=struct('Style','edit','String',60); 
pard.numberOfThetas.position=[5,2];

pard.text5.object=struct('String','overlap in frames','Style','text');
pard.text5.position=[6,1];
pard.numoverlap.object=struct('Style','edit','String',60); 
pard.numoverlap.position=[6,2];
% 
% pard.N0_fit.object=struct('String','N0','Style','radiobutton');
% pard.N0_fit.position=[2,2];
% 
% pard.N0_v.object=struct('String','10','Style','edit');
% pard.N0_v.position=[2,3];
% pard.N0_v.isnumeric=1;
% 
% 
% pard.pmature_fit.object=struct('String','p mature','Style','radiobutton');
% pard.pmature_fit.position=[3,2];
% 
% pard.pmature_v.object=struct('String','.5','Style','edit');
% pard.pmature_v.position=[3,3];
% pard.pmature_v.isnumeric=1;
% 
% 
pard.plugininfo.type='ProcessorPlugin';

end