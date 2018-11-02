classdef intensityTiff<interfaces.SEEvaluationProcessor
    methods
        function obj=intensityTiff(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj,p)
            %  use shift_xy to fit with offset
            if ~isempty(p.layer2_)
            shiftx=p.layer2_.shiftxy_min;
            shifty=p.layer2_.shiftxy_max;
            else
                shiftx=0;
                shifty=0;
            end
            

            pos=obj.site.pos; %position from position of site (center)
            pos(1)=pos(1)-shiftx;
            pos(2)=pos(2)-shifty;
            if p.posfromeval  %use external position instead
                fieldx=p.dxfield;
                fieldy=p.dyfield;
                dx=eval(['obj.site.evaluation.' fieldx]);
                dy=eval(['obj.site.evaluation.' fieldy]);
                pos(1)=pos(1)+dx;
                pos(2)=pos(2)+dy;
               
            end
            
            file=obj.locData.files.file(obj.site.info.filenumber);
            pixcam=file.info.cam_pixelsize_um*1000;
            roi=file.tif(1).info.roi;

            pospixr(1)=round(pos(1)./pixcam(1))-roi(1);
            pospixr(2)=round(pos(2)./pixcam(2))-roi(2);
            
            sfit=round((p.roisize-1)/2);
            
            coim=file.tif(1).image(pospixr(2)-sfit:pospixr(2)+sfit,pospixr(1)-sfit:pospixr(1)+sfit,:);
        
            pr=round(pos([1,2])./pixcam([1,2]));
            rangex=(pr(1)-sfit:pr(1)+sfit)*pixcam(1);
            rangey=(pr(2)-sfit:pr(2)+sfit)*pixcam(2);
                
                
                fixp=[pos(1),pos(2),p.sigmaG,p.mind];
                [X,Y]=meshgrid(rangex,rangey);
                
                
                %find central maximum
                mcentral=max(max(coim(sfit:sfit+2,sfit:sfit+2)));
                
                startim=dgaussforfit([mcentral 0 0 0 0 0],X,Y,[pos(1),pos(2), 250, p.mind]); %single gauss, G1 centered at pos fix, minimum distance p.mind
                
                %determine amplitude and position estimate by subtracting
                %first Gauss.
                [maxrem,ind]=max(coim(:)-startim(:));
                [my,mx]=ind2sub(size(startim),ind);
       
                startp=[mcentral maxrem rangex(mx) rangey(my) 150 0];
                startim=dgaussforfit(startp,X,Y, fixp);
          
                fitp=doublegaussfit(coim,rangex,rangey,startp,fixp);   %fit double gauss, G1 still centered at pos
                
%                 %why still needed? because restricting is only for
%                 drawing, the fitter does not know about it
%                 fo=fitp(3:4);
%                 [fitp(3),fitp(4)]=restrictcoordiantes(fitp(3),fitp(4),pos(1),pos(2),p.mind);
%                 fo-fitp(3:4)
%                 if fitp(3)<pos(1), fitp(3)=min(fitp(3),pos(1)-p.mind); else fitp(3)=max(fitp(3),pos(1)+p.mind); end
%                 if fitp(4)<pos(2), fitp(4)=min(fitp(4),pos(2)-p.mind); else fitp(4)=max(fitp(4),pos(2)+p.mind); end
                
                fitim=dgaussforfit(fitp,X,Y,fixp); %result of the double Gauss fit
                
                %fit poistion of central peak. not used later
%                 startp2=[fitp pos(1) pos(2)];
%                 fixp2=fixp;
%                 fitp2=doublegaussfit2(coim,rangex,rangey,startp2,fixp2);
%                 [fitp2(3),fitp2(4)]=restrictcoordiantes(fitp2(3),fitp2(4),pos(1),pos(2),p.mind);
                
                %fit 3rd Gaussian with fixed xy
                %use residuals to determien start parameter
                dim=coim(:)-fitim(:);
                
                [maxrem,ind]=max(dim);
                [my,mx]=ind2sub(size(fitim),ind);
                
                startp2=[fitp, maxrem,rangex(mx), rangey(my),150];
                fitp3=triplegaussfit(coim,rangex,rangey,startp2,fixp);  %G1 still fixed at center
%                 [fitp3(3),fitp3(4)]=restrictcoordiantes(fitp3(3),fitp3(4),pos(1),pos(2),p.mind);
%                 [fitp3(8),fitp3(9)]=restrictcoordiantes(fitp3(8),fitp3(9),pos(1),pos(2),p.mind);
                
                fitim3=triplegaussforfit(fitp3,X,Y,fixp); %result of 3 Gauss fit
%                 fitimfree=dgaussforfit2(fitp2,X,Y,fixp2);


                % fit again with all 3 Gauss free, but restrict Gauss 1 to
                % be closer than d from pos. 
                fitp3free=triplegaussfitfree(coim,rangex,rangey,[fitp3 pos(1) pos(2)],[p.sigmaG p.mind]);
                
                %chi-squared
                dim=(diff(diff(fitim3-coim,1,1),1,2));
                dvar=var(dim(:))/2;
                
                chi2=sum((fitim3(:)-coim(:)).^2)/dvar;
                chi2red=chi2/(numel(fitim3)-numel(fitp3));
                
                ax=obj.setoutput('image');
                plotfit3(pos,fitp3,ax)
                
                ax=obj.setoutput('freefit');
                plotfit3(fitp3free(11:12),fitp3free(1:10),ax)               
            function plotfit3(pos,fitp3,ax)
                imagesc(rangex,rangey,coim,'Parent',ax)
                 axis(ax,'equal')
                ax.NextPlot='add';
                plot(pos(1),pos(2),'ko','Parent',ax)
                plot(pos(1),pos(2),'k*','Parent',ax)
%                  plot(fitp2(7),fitp2(8),'bx','Parent',ax)
                d=sqrt((fitp3(3)-pos(1))^2+(fitp3(4)-pos(2)).^2);
                sx=(rangex(end)-rangex(1))/2;
                dv=(fitp3(3:4)-pos(1:2));
                mv=max(abs(dv));
                  psh=[mean(rangex) mean(rangey)];
                if mv<=sx
                 plot(fitp3(3),fitp3(4),'k+','Parent',ax)
                 plot(fitp3(3),fitp3(4),'ko','Parent',ax)
                else
                    dv=dv/mv*(sx+0.7*pixcam(1));
                    %%%XXX pixsize conversion: look into
                    
%                     plot(pos(1)+dv(1),pos(2)+dv(2),'k+','Parent',ax)
%                     plot(pos(1)+dv(1),pos(2)+dv(2),'wo','Parent',ax)
                    
                  
                    plot(psh(1)+dv(1),psh(2)+dv(2),'k+','Parent',ax)
                    plot(psh(1)+dv(1),psh(2)+dv(2),'ko','Parent',ax)
                    ax.XLim=[rangex(1)-pixcam(1) rangex(end)+pixcam(1)];
                    ax.YLim=[rangey(1)-pixcam(2) rangey(end)+pixcam(2)];
                end
                
                dv=(fitp3(8:9)-pos(1:2));
                mv=max(abs(dv));
                if mv<=sx
                 plot(fitp3(8),fitp3(9),'kx','Parent',ax)
                 plot(fitp3(8),fitp3(9),'ko','Parent',ax)
                else
                    dv=dv/mv*(sx+0.7*pixcam(1));
                    plot(psh(1)+dv(1),psh(2)+dv(2),'kx','Parent',ax)
                    plot(psh(1)+dv(1),psh(2)+dv(2),'ko','Parent',ax)
                    ax.XLim=[rangex(1)-pixcam(1) rangex(end)+pixcam(1)];
                    ax.YLim=[rangey(1)-pixcam(2) rangey(end)+pixcam(2)];
                end
                
                ttxt=['*: A1=' num2str(fitp3(1),'%5.0f') ', +: A2=' num2str(fitp3(2),'%5.0f') ', x: A3=' num2str(fitp3(7),'%5.0f') ', d=' num2str(d,'%5.0f') ', chi2=' num2str(chi2red,'%5.2f')];
                title(ttxt,'Parent',ax)
                axis(ax,'equal')
                ax.NextPlot='replace';
            end
                ax2=obj.setoutput('fit');
                
                imagesc(rangex,rangey,fitim3,'Parent',ax2)
                title(fitp(1:2),'Parent',ax2)
                axis(ax2,'equal')
                ax3=obj.setoutput('residuals');
                
                imagesc(rangex,rangey,coim-fitim3,'Parent',ax3) 
                axis(ax3,'equal')
                
                ax4=obj.setoutput('startim');
                
                
                imagesc(rangex,rangey,startim,'Parent',ax4)   
                axis(ax4,'equal')
                
                
                out.G3Amplitude1=fitp3(1);
                out.G3Amplitude2=fitp3(2);
                out.G3Amplitude3=fitp3(7);
                
                out.G3freeAmplitude1=fitp3free(1);
                out.G3freeAmplitude2=fitp3free(2);
                out.G3freeAmplitude3=fitp3free(7);
                
                out.G3freexcent=fitp3free(11);
                out.G3freeycent=fitp3free(12);
                
%                 out.Amplitudefreexy1=fitp2(1);
%                 out.Amplitudefreexy2=fitp2(2);
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end

function pard=guidef
pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi','layer2_'};
pard.t1.object=struct('Style','text','String','ROI (pixels)');
pard.t1.position=[1,1];
pard.t1.Width=2;

pard.roisize.object=struct('Style','edit','String','7');
pard.roisize.position=[1,4];

pard.t2.object=struct('Style','text','String','minimal distance 2nd Gauss (nm)');
pard.t2.position=[2,1];
pard.t2.Width=3;

pard.mind.object=struct('Style','edit','String','350');
pard.mind.position=[2,4];

pard.t3.object=struct('Style','text','String','sigma first Gauss (nm)');
pard.t3.position=[3,1];
pard.t3.Width=2;

pard.sigmaG.object=struct('Style','edit','String','150');
pard.sigmaG.position=[3,4];
pard.plugininfo.type='ROI_Evaluate';

pard.posfromeval.object=struct('Style','checkbox','String','use position from evaluation field: site.evaluation.');
pard.posfromeval.position=[5,1];
pard.posfromeval.Width=4;
pard.t4.object=struct('Style','text','String','dx0');
pard.t4.position=[6,1];
pard.t4.Width=.5;
pard.t5.object=struct('Style','text','String','dy0');
pard.t5.position=[7,1];
pard.t5.Width=.5;

pard.dxfield.object=struct('Style','edit','String','CME2DRing.imfit.x0');
pard.dxfield.position=[6,1.5];
pard.dxfield.Width=3.5;
pard.dyfield.object=struct('Style','edit','String','CME2DRing.imfit.y0');
pard.dyfield.position=[7,1.5];
pard.dyfield.Width=3.5;

% 
% pard.fitxy.object=struct('Style','checkbox','String','Fit also position of DL peaks');
% pard.fitxy.position=[6,1];
% pard.fitxy.Width=4;

pard.plugininfo.type='ROI_Evaluate';

end



function fit=doublegaussfit(img,rangex,rangey,startp,fixp)
[X,Y]=meshgrid(double(rangex),double(rangey));

dx=rangex(end)-rangex(1);
meanposx=mean(rangex); meanposy=mean(rangey);
n=1;
lb=double([0 0 meanposx-n*dx meanposy-n*dx 150 0]);ub=double([inf inf meanposx+n*dx meanposy+n*dx  750 inf]);
fit=lsqnonlin(@dgaussforfiterr,double(startp),lb,ub,[],double(img),X,Y,fixp);

[fit(3),fit(4)]=restrictcoordiantes(fit(3),fit(4),fixp(1),fixp(2),fixp(4));
end
function out=dgaussforfiterr(fitp,img,X,Y,fixp)
out=dgaussforfit(fitp,X,Y,fixp)-img;
end
function out=dgaussforfit(fitp,X,Y,fixp)
% fitp: A1 A2 x2 y2 sigma2 bg
%fixp x1, y1, sigma1

A1=fitp(1); A2=fitp(2); x2=fitp(3); y2=fitp(4); sigma2=fitp(5); bg=fitp(6);

d=fixp(4);
x1=fixp(1);y1=fixp(2); sigma1=fixp(3);
%adjsut x1, y1
% if x2<x1, x2=min(x2,x1-d); else x2=max(x2,x1+d); end
% if y2<y1, y2=min(y2,y1-d); else y2=max(y2,y1+d); end
[x2,y2]=restrictcoordiantes(x2,y2,x1,y1,d);
out=double(A1* exp( - ((X-x1).^2+(Y-y1).^2)/sigma1^2/2)+A2* exp( - ((X-x2).^2+(Y-y2).^2)/sigma2^2/2)+bg);

end

function fit=triplegaussfit(img,rangex,rangey,startp,fixp)
[X,Y]=meshgrid(double(rangex),double(rangey));

dx=rangex(end)-rangex(1);
meanposx=mean(rangex); meanposy=mean(rangey);

% change n to determine how far 2nd and 3rd Gaussians can move outside, n
% is maximum distance from center in [ROI size]s
n=.5;
lb=double([0 0 meanposx-n*dx meanposy-n*dx 150 0 0 meanposx-n*dx meanposy-n*dx  150]);ub=double([inf inf meanposx+n*dx meanposy+n*dx  250 inf inf meanposx+n*dx meanposy+n*dx 250]);
startp(startp>ub)=ub(startp>ub);
startp(startp<lb)=lb(startp<lb);
fit=lsqnonlin(@triplegaussforfiterr,double(startp),lb,ub,[],double(img),X,Y,fixp);
[fit(3),fit(4)]=restrictcoordiantes(fit(3),fit(4),fixp(1),fixp(2),fixp(4));
[fit(8),fit(9)]=restrictcoordiantes(fit(8),fit(9),fixp(1),fixp(2),fixp(4));
end
function out=triplegaussforfiterr(fitp,img,X,Y,fixp)
out=triplegaussforfit(fitp,X,Y,fixp)-img;
end
function out=triplegaussforfit(fitp,X,Y,fixp)
% fitp: A1 A2 x2 y2 sigma2 bg A3 x3 y3 sigma3
%fixp x1, y1, sigma1

% A1=fitp(1); A2=fitp(2); x2=fitp(3); y2=fitp(4); sigma2=fitp(5); bg=fitp(6);

A3=fitp(7); x3=fitp(8); y3=fitp(9); sigma3=fitp(10);

d=fixp(4);
x1=fixp(1);y1=fixp(2); 
% igma1=fixp(3);
%adjsut x1, y1
% if x2<x1, x2=min(x2,x1-d); else x2=max(x2,x1+d); end
% if y2<y1, y2=min(y2,y1-d); else y2=max(y2,y1+d); end
% [x2,y2]=restrictcoordiantes(x2,y2,x1,y1,d);
[x3,y3]=restrictcoordiantes(x3,y3,x1,y1,d);

out=dgaussforfit(fitp(1:6),X,Y,fixp)+double(A3* exp( - ((X-x3).^2+(Y-y3).^2)/sigma3^2/2));
% out=double(A1* exp( - ((X-x1).^2+(Y-y1).^2)/sigma1^2/2)+A2* exp( - ((X-x2).^2+(Y-y2).^2)/sigma2^2/2)+bg);

end


function out=triplegaussforfitfree(fitp,X,Y,fixp)
% fitp: A1 A2 x2 y2 sigma2 bg A3 x3 y3 sigma3
%fixp x1, y1, sigma1

A1=fitp(1); A2=fitp(2); x2=fitp(3); y2=fitp(4); sigma2=fitp(5); bg=fitp(6);

A3=fitp(7); x3=fitp(8); y3=fitp(9); sigma3=fitp(10);

d=fixp(1);
% d=fixp(4);

x1=fitp(11);y1=fitp(12); 

fixp2gauss=[x1,y1,fixp(2), d];
% igma1=fixp(3);
%adjsut x1, y1
% if x2<x1, x2=min(x2,x1-d); else x2=max(x2,x1+d); end
% if y2<y1, y2=min(y2,y1-d); else y2=max(y2,y1+d); end
% [x2,y2]=restrictcoordiantes(x2,y2,x1,y1,d);
[x3,y3]=restrictcoordiantes(x3,y3,x1,y1,d);
% [x2,y2]=restrictcoordiantes(x2,y2,x1,y1,d);

out=dgaussforfit(fitp(1:6),X,Y,fixp2gauss)+double(A3* exp( - ((X-x3).^2+(Y-y3).^2)/sigma3^2/2));
% out=double(A1* exp( - ((X-x1).^2+(Y-y1).^2)/sigma1^2/2)+A2* exp( - ((X-x2).^2+(Y-y2).^2)/sigma2^2/2)+bg);

end

function out=triplegaussforfiterrfree(fitp,img,X,Y,fixp)
out=triplegaussforfitfree(fitp,X,Y,fixp)-img;
end

function fit=triplegaussfitfree(img,rangex,rangey,startp,fixp)
d=fixp(2);
[X,Y]=meshgrid(double(rangex),double(rangey));

dx=rangex(end)-rangex(1);
meanposx=mean(rangex); meanposy=mean(rangey);

% change n to determine how far 2nd and 3rd Gaussians can move outside, n
% is maximum distance from center in [ROI size]s
n=.5;
lb=double([0 0 meanposx-n*dx meanposy-n*dx 150 0 0 meanposx-n*dx meanposy-n*dx  150 meanposx-d meanposy-d]);ub=double([inf inf meanposx+n*dx meanposy+n*dx  250 inf inf meanposx+n*dx meanposy+n*dx 250 meanposx+d meanposy+d]);
fit=lsqnonlin(@triplegaussforfiterrfree,double(startp),lb,ub,[],double(img),X,Y,fixp);
[fit(3),fit(4)]=restrictcoordiantes(fit(3),fit(4),fit(11),fit(12),fixp(2));
[fit(8),fit(9)]=restrictcoordiantes(fit(8),fit(9),fit(11),fit(12),fixp(2));
end
% function fit=doublegaussfit2(img,rangex,rangey,startp,fixp)
% [X,Y]=meshgrid(double(rangex),double(rangey));
% dx=(rangex(2)-rangex(1))/4;
% lb=double([0 0 -inf -inf 150 0 rangex(1)+dx rangey(1)+dx]);ub=double([inf inf inf inf 750 inf rangex(end)-dx rangey(end)-dx]);
% fit=lsqnonlin(@dgaussforfiterr2,double(startp),lb,ub,[],double(img),X,Y,fixp);
% 
% end
% function out=dgaussforfiterr2(fitp,img,X,Y,fixp)
% out=dgaussforfit2(fitp,X,Y,fixp)-img;
% end
% function out=dgaussforfit2(fitp,X,Y,fixp)
% % fitp: A1 A2 x2 y2 sigma2 bg
% %fixp x1, y1, sigma1
% 
% A1=fitp(1); A2=fitp(2); x2=fitp(3); y2=fitp(4); sigma2=fitp(5); bg=fitp(6);
% 
% d=fixp(4);
% x1=fitp(7);y1=fitp(8); sigma1=fixp(3);
% %adjsut x1, y1
% % if x2<x1, x2=min(x2,x1-d); else x2=max(x2,x1+d); end
% % if y2<y1, y2=min(y2,y1-d); else y2=max(y2,y1+d); end
% [x2,y2]=restrictcoordiantes(x2,y2,x1,y1,d);
% out=double(A1* exp( - ((X-x1).^2+(Y-y1).^2)/sigma1^2/2)+A2* exp( - ((X-x2).^2+(Y-y2).^2)/sigma2^2/2)+bg);
% 
% end

function [x2,y2]=restrictcoordiantes(x2,y2,x1,y1,d)
dvec=[x2-x1;y2-y1];
dist=norm(dvec);
% dist=sqrt((x2-x1)^2+(y2-y1)^2);
if dist<d
    xnew=[x1;y1]+dvec/dist*d;
    x2=xnew(1);
    y2=xnew(2);

end
end