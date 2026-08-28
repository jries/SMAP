classdef correlationAnalysis<interfaces.SEEvaluationProcessor
%     Calcualtes the shift between layer 1 and layer 2 based on
%     cross-correlation.
    properties

    end
    methods
        function obj=correlationAnalysis(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj,p)
%             obj.line=p.lineselect.selection;
           
%             if isempty(obj.site.image.layers)
%                 warndlg('check "keep temp images" in the Evaluation GUI');
%                 out=[];
%                 return
%             end
            
%             maxprecision=0.1;% nm
            maxprecision=1;% nm

            %get coordinates
            layers=find(p.sr_layerson);
            locs=obj.getloc({'xnm','ynm','znm','locprecnm','layer'},'layer',layers,'size','freeroi');
            %rotate with respect to line
            
            lineangle=obj.site.annotation.(p.lineselect.selection).angle;
            linepos=obj.site.annotation.(p.lineselect.selection).pos;
            if all(linepos==0)
                out=[];
                return
            end
            linelen=sqrt(sum((linepos(2,:)-linepos(1,:)).^2))*1000;
            mpos=mean(linepos,1)*1000;
            
            [xr,yr]=rotcoorddeg(locs.xnm-mpos(1),locs.ynm-mpos(2),lineangle);
            inroixy=abs(yr)<=p.width/2 & abs(xr)<=linelen/2;
            inz=locs.znm>p.zrange(1) & locs.znm<p.zrange(2);
            zdist=50;
            zpos=p.zrange(1):zdist:p.zrange(2);
            zwidth=100;
            ydist=50;
            ywidth=100;
            ypos=-p.width/2:ydist:p.width/2;
            [acimxy,nx,ny,imxy]=ac2dslices(xr(inroixy&inz),yr(inroixy&inz),locs.znm(inroixy&inz),[-linelen/2 linelen/2], [-1 1]*p.width/2, p.pixrec, zpos, zwidth);
            [acimxz,nxz,nz,imxz]=ac2dslices(xr(inroixy&inz),locs.znm(inroixy&inz),yr(inroixy&inz),[-linelen/2 linelen/2], p.zrange, p.pixrec, ypos, ywidth);
            acimxy(end,:)=max(acimxy(:));
            % axac=obj.setoutput("cc2d");
            % sim=size(acimxz,1);
            acimt=vertcat(acimxy,acimxz);acimt=acimt/max(acimt(:));
            % imagesc(axac,vertcat(acimxy,acimxz(round(sim*.25):round(sim*.75),:)))
            % imagesc(axac,vertcat(acimxy,acimxz))
            % axis(axac,"equal")
            imxy=imgaussfilt(imxy,1);imxz=imgaussfilt(imxz,01);
            imxy=imxy/quantile(imxy(:),.99);
            imxy(end,:)=1;
            imxz=imxz/quantile(imxz(:),.99);
            imt=vertcat(imxy,imxz);
            imt(imt>1)=1;
            axim=obj.setoutput("img");
            imt(:,end)=max(imt(:));

            imdisp=horzcat(imt,acimt);
            nxd=(1:size(imdisp,2))*p.pixrec;
            nyd=(1:size(imdisp,1))*p.pixrec;

            scalebarlen=round(100/p.pixrec);
            imdisp(end-2,4:4+scalebarlen-1)=max(imdisp(:));
            imagesc(axim,nxd,nyd,imdisp)
            axis(axim,"equal")
            % colormap(axim,"hot")
            
            wins=200/p.pixrec/2;
            midpxy=ceil(size(acimxy)/2);
            acxyzoom=acimxy(midpxy(1)-wins:midpxy(1)+wins,midpxy(2)-wins:midpxy(2)+wins);
            midpxz=ceil(size(acimxz)/2);
            acxzzoom=acimxz(midpxz(1)-wins:midpxz(1)+wins,midpxz(2)-wins:midpxz(2)+wins);
            acxyzoom(:,end)=max(acxyzoom(:));


            axac=obj.setoutput("cc2dz");

            imdispz=horzcat(acxyzoom,acxzzoom);
            nxd=(1:size(imdispz,2))*p.pixrec;
            nyd=(1:size(imdispz,1))*p.pixrec;

            scalebarlen=round(100/p.pixrec);
            imdispz(end-2,4:4+scalebarlen-1)=max(imdispz(:));

            imagesc(axac,nxd,nyd,imdispz)
            axis(axac,"equal")


            axprofiles=obj.setoutput("profiles");
            
            winprof=round(5/2);
            midp=ceil(size(acimxy)/2);        
            px=mean(acimxy(midp(1)-winprof:midp(1)+winprof,:),1);
            midpz=ceil(size(acimxz)/2);        
            pz=mean(acimxz(midpz(1)-winprof:midpz(1)+winprof,:),1);

            % py=mean(acimxy,2);
            % pz=mean(acimxz,2);
            % pzx=mean(acimxz,1);
            hold(axprofiles,"off")
            plot(axprofiles,nx(ceil(end/2):end),px(ceil(end/2):end));%,ny(ceil(end/2):end),py(ceil(end/2):end),nz(ceil(end/2):end),pz(ceil(end/2):end),nx(ceil(end/2):end),pzx(ceil(end/2):end))
             hold(axprofiles,"on")
            plot(axprofiles,nx(ceil(end/2):end),pz(ceil(end/2):end));
            % legend(axprofiles,"x","y","z","x(z)")

            out.cc.acxy=acimxy;
            out.cc.acxz=acimxz;
        end
        function pard=guidef(obj)
            pard=guidef;
        end

    end
end


function [acim,nx,ny,imtot]=ac2dslices(x,y,z,xrange, yrange, pixelsize, zpos, zwidth)
    nx=xrange(1):pixelsize:xrange(2);
    ny=yrange(1):pixelsize:yrange(2);
    acim=zeros(length(ny)-1,length(nx)-1);
    imtot=zeros(length(ny)-1,length(nx)-1);
    for k=1:length(zpos)
        inz=z>=zpos(k) & z<zpos(k)+zwidth;
        img=histcounts2(y(inz),x(inz),ny,nx);
        acim=acim+accorr2fft(img);
        imtot=imtot+img;
    end
    acim=acim/length(zpos);
    acim(ceil(length(ny)/2),ceil(length(nx)/2))=NaN;        
    nx=nx(1:end-1);ny=ny(1:end-1);
end

function out=accorr2fft(in1)
    fi1=fft2(in1); 
    cc=fi1.*conj(fi1);
    out=fftshift(ifft2(cc));
end


% function [yshift, yProfile, gaussFit, bg] = auto_corr_analysis(sxf, p, maxprecision)
%     % Resize to one pixel per shift
%     [~,linind] = max(sxf(:));
%     [xm,ym]=ind2sub(size(sxf),linind);
%     yProfile = sxf(xm,:);
%     [~,linind] = max(yProfile);
%     yProfile(linind) = yProfile(linind+1);
% 
%     yProfile=imresize(yProfile,[1 length(yProfile)*p.pixrec/maxprecision],'cubic');
% %     sxhr = sxf;
% %     maxprecision = 1;
% %     [~,linind]=max(sxhr(:));
% %     [xm,ym]=ind2sub(size(sxhr),linind);
% % out.dxline=(xm-size(sx12hr,1)/2)*maxprecision;
% %     yProfile = sxhr(xm,:);
%     yProfile = yProfile-min(yProfile); 
% %             yMax_val = max(yProfile);
% %     gauss2 = 'a1*exp(-((x-b)/c1)^2) + a2*exp(-((x-b)/c2)^2)+d';
%     gauss2 = 'a1*exp(-((x-b)/c1)^2) + a2*exp(-((x-b)/c2)^2)';
%     yshift = (-length(yProfile)/2:length(yProfile)/2).*maxprecision;
%     yshift = yshift(1:end-1);
% %     gaussFit = fit(yshift',yProfile', gauss2, 'StartPoint',[max(yProfile)*0.9 max(yProfile)*0.1 0 10 100 0]);
%     gaussFit = fit(yshift',yProfile', gauss2, 'StartPoint',[max(yProfile)*0.9 max(yProfile)*0.1 0 25 100]);
% %     bg = gaussFit.a2+gaussFit.d;
%     bg = gaussFit.a2;
% end
% 
% function [yshift, yProfile, gaussFit, bg] = auto_corr_analysis(sxf, p, maxprecision)
%     sxhr=imresize(sxf,p.pixrec/maxprecision,'cubic');
%     % Resize to one pixel per shift
% %     maxprecision = 1;
%     [~,linind]=max(sxhr(:));
%     [xm,ym]=ind2sub(size(sxhr),linind);
% % out.dxline=(xm-size(sx12hr,1)/2)*maxprecision;
%     yProfile = sxhr(xm,:);
%     yProfile = yProfile-min(yProfile); 
% %             yMax_val = max(yProfile);
% %     gauss2 = 'a1*exp(-((x-b)/c1)^2) + a2*exp(-((x-b)/c2)^2)+d';
%     gauss2 = 'a1*exp(-((x-b)/c1)^2) + a2*exp(-((x-b)/c2)^2)';
%     yshift = (-size(sxhr,1)/2:size(sxhr,1)/2).*maxprecision;
%     yshift = yshift(1:end-1);
% %     gaussFit = fit(yshift',yProfile', gauss2, 'StartPoint',[max(yProfile)*0.9 max(yProfile)*0.1 0 10 100 0]);
%     gaussFit = fit(yshift',yProfile', gauss2, 'StartPoint',[max(yProfile)*0.9 max(yProfile)*0.1 0 25 100]);
% %     bg = gaussFit.a2+gaussFit.d;
%     bg = gaussFit.a2;
% end
function pard=guidef
pard.inputParameters={'se_sitepixelsize','numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi','layer2_'};

pard.lineselect.object=struct('Style','popupmenu','String',{{'line1','line2','rotationpos'}},'Value',3);
pard.lineselect.position=[1,1];
pard.lineselect.Width=4;

%maximum shift
% fit window

pard.widtht.object=struct('Style','text','String','line width (nm)');
pard.widtht.position=[2,1];
pard.widtht.Width=1;

pard.width.object=struct('Style','edit','String','250');
pard.width.position=[2,2];
pard.width.Width=.5;

pard.zranget.object=struct('Style','text','String','z range (nm)');
pard.zranget.position=[2,3];
pard.zranget.Width=1;

pard.zrange.object=struct('Style','edit','String','-400 400');
pard.zrange.position=[2,4];
pard.zrange.Width=1;

pard.pixrect.object=struct('Style','text','String','pixel size (nm)');
pard.pixrect.position=[3,1];
pard.pixrect.Width=1;

pard.pixrec.object=struct('Style','edit','String','5');
pard.pixrec.position=[3,2];
pard.pixrec.Width=.5;

% 
% pard.shiftval.object=struct('Style','edit','String','0');
% pard.shiftval.position=[2,3];
% pard.shiftval.Width=2;

pard.plugininfo.type='ROI_Evaluate';
pard.plugininfo.description='Calcualtes the shift between layer 1 and layer 2 based on cross-correlation.';
end


