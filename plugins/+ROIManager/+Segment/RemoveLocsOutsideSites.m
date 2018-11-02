classdef RemoveLocsOutsideSites<interfaces.DialogProcessor&interfaces.SEProcessor
    methods
        function obj=RemoveLocsOutsideSites(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'se_sitefov','se_cellpixelsize','se_siteroi'};
        end
        
        function out=run(obj,p)  
          removeLocsSites(obj,p)
          out=[];
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end




function pard=guidef
% pard.t1.object=struct('String','cutoffmask','Style','text');
% pard.t1.position=[2,1];
% pard.cutoffmask.object=struct('String','5','Style','edit');
% pard.cutoffmask.position=[2,2];
% 
% pard.t2.object=struct('String','filtersize (nm)','Style','text');
% pard.t2.position=[1,1];
% pard.sigma1.object=struct('String','200','Style','edit');
% pard.sigma1.position=[1,2];
% 
% pard.t3.object=struct('String','rim / radius','Style','text');
% pard.t3.position=[4,1];
% pard.rim.object=struct('String','0.2','Style','edit');
% pard.rim.position=[4,2];
% 
% pard.t9.object=struct('String','iterations rim dilation','Style','text');
% pard.t9.position=[3,1];
% pard.iterations.object=struct('String','3','Style','edit');
% pard.iterations.position=[3,2];
% 
% 
% pard.t4.object=struct('String','cutoff maxima','Style','text');
% pard.t4.position=[2,3];
% pard.cutoff.object=struct('String','20','Style','edit');
% pard.cutoff.position=[2,4];
% 
% pard.t5.object=struct('String','max pos corr (nm)','Style','text');
% pard.t5.position=[3,3];
% pard.maxjump.object=struct('String','75','Style','edit');
% pard.maxjump.position=[3,4];
% 
% pard.t6.object=struct('String','size range (std) (nm)','Style','text');
% pard.t6.position=[4,3];
% pard.sizerange.object=struct('String','20 200','Style','edit');
% pard.sizerange.position=[4,4];
% 
% pard.t7.object=struct('String','filter findmax (nm)','Style','text');
% pard.t7.position=[1,3];
% pard.sigmafind.object=struct('String','100','Style','edit');
% pard.sigmafind.position=[1,4];
% 
% pard.preview.object=struct('String','preview','Style','checkbox','Value',1);
% pard.preview.position=[6,1];
% 
% pard.t11.object=struct('String','Works best with: PSF filtered, locprec filter larger than usual (~30-60 nm), grouped','Style','text');
% pard.t11.position=[8,1];
% pard.t11.Width=4;
pard.plugininfo.type='ROI_Analyze';
end

function removeLocsSites(obj,p)
% global se
se=obj.SE;
locData=obj.locData;
indin=false(size(locData.loc.xnm));
sites=se.sites;
for k=1:length(sites)
    site=sites(k);
    indin=indin|(locData.loc.filenumber==site.info.filenumber & mywithin(locData.loc.xnm,[site.pos(1)-p.se_sitefov/2,p.se_sitefov],locData.loc.ynm,[site.pos(2)-p.se_sitefov/2,p.se_sitefov]));
end

obj.locData.removelocs(~indin);
obj.locData.regroup;


% cells=se.cells;
% pr=p.se_cellpixelsize;
% sigma1=p.sigma1/pr;
% sigmafind=p.sigmafind/pr;
% % sigma1=25;
% 
% rimrel=p.rim;
% % rim=20;
% if length(p.sizerange)==1;
%     p.sizerange(2)=inf;
% end
% hgauss=fspecial('gaussian',ceil(sigma1*5),sigma1);
% % hgauss=fspecial('disk',sigma1);
% hgmax=max(hgauss(:));
% hgauss3=fspecial('gaussian',15,2);
% hgauss2=fspecial('gaussian',ceil(sigmafind*5),sigmafind);
% cutoff=p.cutoff;
% cutoffmax=p.cutoffmask;
% 
% if p.preview
%     rangec=se.indexFromID(se.cells,se.currentcell.ID);
% else
%     rangec=1:length(cells);
% end
% 
% struce = strel('disk',5);
% struce2 = strel('disk',1);
% 
%     hf=initaxis(p.resultstabgroup,'segmentation');
%     
%     
% for kc=rangec
%     cell=se.cells(kc);
%     locs=locData.getloc({'xnm','ynm'},'layer',1,'position',[cell.pos(1:2) p.se_cellfov p.se_cellfov],'filenumber',cell.info.filenumber);
%     rangex=[cell.pos(1)-p.se_cellfov/2 , cell.pos(1)+p.se_cellfov/2 ];
%     rangey=[cell.pos(2)-p.se_cellfov/2 , cell.pos(2)+p.se_cellfov/2 ];
%     imh=(myhist2(locs.xnm,locs.ynm,pr,pr, rangex,rangey))';
%     imfs=imfilter(sqrt(imh),hgauss);
%     imf=imfilter(imh,hgauss3);
%     imbw=imfs>hgmax*cutoffmax;
%     imbw=imclose(bwareafilt(imbw,1),struce);
%     
%     xrel=(locs.xnm-rangex(1))/pr;
%     yrel=(locs.ynm-rangey(1))/pr;
% %     indin=withinmask(imbw,yrel,xrel);
%     che=imbw;
%     for k=1:p.iterations
%         indin=withinmask(che,yrel,xrel);
%         ims=(myhist2(locs.xnm(indin),locs.ynm(indin),pr,pr, rangex,rangey))';
%         ch=bwconvhull(ims>0);
%         che=imerode(ch,struce2);
%     end
% %     inm=sum(indin(:))/length(indin);
%     %rim relative 
%     areach=sum(ch(:));
%     radiusch=sqrt(areach/pi);
%     rim=radiusch*rimrel;
%     struce3= strel('disk',round(rim));
%     che=imerode(ch,struce3);
%     imh2=imh;imh2(~che)=0;
%     imincell=imfilter(imh2,hgauss2);
%      maxima=maximumfindcall(imincell);
%      
%       indco=maxima(:,3)>cutoff*max(hgauss2(:));
%     
%     cx=maxima(indco,1)*pr+rangex(1);
%     cy=maxima(indco,2)*pr+rangey(1);
%      
%     s=size(imh);
% %     implot=zeros(s(1),s(2),3);
%     map=hot(256);
%     map2=gray(256);
%     norm=myquantilefast(imf,[.98,.9999]);
%     imbg=imf/norm(1);imbg(imbg>1)=1;
%     implot=ind2rgb(uint8(imf/norm(2)*255),map);%+0*ind2rgb(uint8(imf/norm(1)*255)/4,map2);
%     implot(:,:,2)=implot(:,:,2)+bwperim((imbw))+imbg/4;
% %     implot(:,:,2)=imf/max(imf(:))*2;
%     implot(:,:,3)=implot(:,:,3)+bwperim(che);
%     implot(:,:,1)=implot(:,:,1)+bwperim(ch);
%     
%     hold off
%     imagesc(rangex,rangey,implot,'Parent',hf)
%       hold on
%     plot(cx,cy,'ko','Parent',hf)
%     title(hf,kc);
%     
%     ind=1;
%     ckeepx=[];ckeepy=[];
%     for ks=1:length(cx)
%         %calculate average position
%         pos=[cx(ks),cy(ks)];
%         locs=locData.getloc({'xnm','ynm'},'layer',1,'position',[pos(1:2) p.se_siteroi p.se_siteroi],'filenumber',cell.info.filenumber);
%         posnew=[mean(locs.xnm), mean(locs.ynm)];
%         
%         stdposx=std(locs.xnm);
%         stdposy=std(locs.ynm);
%         inpos=stdposx>p.sizerange(1)&&stdposx<p.sizerange(2);
%         inpos=inpos&(stdposy>p.sizerange(1)&&stdposy<p.sizerange(2));
%         
%         
%         if sum((pos-posnew).^2)<p.maxjump^2  &&inpos
%         %if htat is too far from current position: was outside rim: remove
%             ckeepx(ind)=posnew(1);ckeepy(ind)=posnew(2);ind=ind+1;
%         
%             currentsite=interfaces.SEsites;
%             currentsite.pos=posnew;
%             currentsite.ID=0;
% 
%         %     currentcell.sePar=obj.SE.sePar;
%             currentsite.info.cell=cell.ID;
%             currentsite.info.filenumber=cell.info.filenumber;
%             if ~p.preview
%                 obj.SE.addSite(currentsite);
%             end
%         end
%     end
%     
%         plot(ckeepx,ckeepy,'kx','Parent',hf)
%         drawnow
%         
%     
% end
% 
% se.processors.preview.updateSitelist
end
