classdef shiftCrossCorrelation<interfaces.SEEvaluationProcessor
%     Calcualtes the shift between layer 1 and layer 2 based on
%     cross-correlation.
    properties

    end
    methods
        function obj=shiftCrossCorrelation(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj,p)
%             obj.line=p.lineselect.selection;
           
%             if isempty(obj.site.image.layers)
%                 warndlg('check "keep temp images" in the Evaluation GUI');
%                 out=[];
%                 return
%             end
            
            maxprecision=0.1;% nm
       
            %get coordinates
            layers=find(p.sr_layerson);
            locs1=obj.getloc({'xnm','ynm','locprecnm','layer'},'layer',layers(1),'size','freeroi');
            locs2=obj.getloc({'xnm','ynm','locprecnm','layer'},'layer',layers(2),'size','freeroi');
            %rotate with respect to line
            
            lineangle=obj.site.annotation.(p.lineselect.selection).angle;
            
            [xr1,yr1]=rotcoorddeg(locs1.xnm-obj.site.pos(1),locs1.ynm-obj.site.pos(2),lineangle);
            [xr2,yr2]=rotcoorddeg(locs2.xnm-obj.site.pos(1),locs2.ynm-obj.site.pos(2),lineangle);
            n=-p.se_siteroi/2:p.pixrec:p.se_siteroi/2;
            img1=histcounts2(xr1,yr1,n,n);
            img2=histcounts2(xr2,yr2,n,n);
            satpix=10;
            qp=1-satpix/size(img1,1)/size(img1,2);
            qp=max(qp,0.99);
            q1=max(2,ceil(quantile(img1(:),qp)));
            img1(img1>q1)=q1;
            q2=max(2,ceil(quantile(img2(:),qp)));
            img2(img2>q2)=q2;           
            %saturate to reduce impact of bright spots
            
%             img1=img1-mean(img1(:));
%             img2=img2-mean(img2(:));

            if obj.display
                ax1=obj.setoutput('images');
                imagesc(ax1,horzcat(img1,img2));
            end
            %make mask
            sigma=p.fitsig/p.pixrec;
            
            h=fspecial('gaussian',ceil(3*sigma),sigma);
            x1=xcorr2fft(img1,img1);
            x2=xcorr2fft(img2,img2);
            x12=xcorr2fft(img1,img2);
            x12f=filter2(h,x12);
            
            [v,I] = sort(x1(:));
            x1rmp = x1;
            x1rmp(I(1)) = v(2);
            x1f=filter2(h,x1rmp);

            [v,I] = sort(x2(:));
            x2rmp = x2;
            x2rmp(I(1)) = v(2);
            x2f=filter2(h,x2rmp);
             maxshiftwin=ceil(p.maxshift/p.pixrec)+ceil((3*sigma-1)/2);
             midp=floor(size(x12)/2+1);
             sx1=x1(midp(1)-maxshiftwin:midp(1)+maxshiftwin,midp(2)-maxshiftwin:midp(2)+maxshiftwin);
             sx2=x2(midp(1)-maxshiftwin:midp(1)+maxshiftwin,midp(2)-maxshiftwin:midp(2)+maxshiftwin);
             sx12=x12(midp(1)-maxshiftwin:midp(1)+maxshiftwin,midp(2)-maxshiftwin:midp(2)+maxshiftwin);
             sx12f=x12f(midp(1)-maxshiftwin:midp(1)+maxshiftwin,midp(2)-maxshiftwin:midp(2)+maxshiftwin);
             sx1f=x1f(midp(1)-maxshiftwin:midp(1)+maxshiftwin,midp(2)-maxshiftwin:midp(2)+maxshiftwin);
             sx2f=x2f(midp(1)-maxshiftwin:midp(1)+maxshiftwin,midp(2)-maxshiftwin:midp(2)+maxshiftwin);
             
             sx12hr=imresize(sx12f,p.pixrec/maxprecision,'cubic');
             [~,linind]=max(sx12hr(:));
             [xm,ym]=ind2sub(size(sx12hr),linind);
             if obj.display
                 ax2=obj.setoutput('corr');
                 imagesc(ax2,[sx1/max(sx1(:)),sx2/max(sx2(:));sx12/max(sx12(:)),sx12f/max(sx12f(:))]);                  
             end

             %cut out center and fit to get precise center. First output
             % project along line, summ up perpendicular direction -> cc
             % along line for further averaging.
            out.xcorr=sx12;
            out.dxline=(xm-size(sx12hr,1)/2)*maxprecision;
            out.dyline=(ym-size(sx12hr,2)/2)*maxprecision;

            out.xr1 = xr1;
            out.yr1 = yr1;
            out.xr2 = xr2;
            out.yr2 = yr2;
            

            % YW added:
            % get the FWHM of the peak of the xcross perpendicular to the
            % kinetochore axis.
%             yProfile = sx12hr(xm,:);
%             yProfile = yProfile-min(yProfile); 
%             yMax_val = max(yProfile);
%             i1=find(yProfile>yMax_val/2,1,'first');
%             i2=find(yProfile>yMax_val/2,1,'last');
%             fwhm = (i2-i1)*maxprecision;

%             ax3=obj.setoutput('fwhm');
%             yshift = (-size(sx12hr,1)/2:size(sx12hr,1)/2)*maxprecision;
%             plot(ax3, yshift(1:end-1)+0.5, yProfile,'r');
%             hold(ax3, 'on')
%             yline(ax3, yMax_val/2)
%             title(ax3, ['FWHM = ' num2str(fwhm) ' nm'])
%             xlabel(ax3, 'y shift (nm)')
%             ylabel(ax3, 'xcorr')
%             hold(ax3, 'off')

%             out.fwhm12 = fwhm;

            [yshift1, yProfile1, gaussFit1, bg1] = auto_corr_analysis(sx1f, p, maxprecision);
            [yshift2, yProfile2, gaussFit2, bg2] = auto_corr_analysis(sx2f, p, maxprecision);
            
            if obj.display
                val = feval(gaussFit1, -p.maxshift:p.maxshift);
    %             i1=find(yProfile>yMax_val/2,1,'first');
    %             i2=find(yProfile>yMax_val/2,1,'last');
    %             fwhm = (i2-i1)*maxprecision;
    %             size_oneK = fwhm - 2.355*mean(locs1.locprecnm);
                ax4=obj.setoutput('fwhm_auto_l1');
                plot(ax4, yshift1, yProfile1,'r');
                hold(ax4, 'on')
    %             yline(ax4, yMax_val/2)
                plot(ax4, -p.maxshift:p.maxshift, val);
                yline(ax4, bg1);
    %             title(ax4, ['FWHM = ' num2str(fwhm) ' nm; width_k = ' num2str(size_oneK)])
                xlabel(ax4, 'y shift (nm)')
                ylabel(ax4, 'auto corr')
                hold(ax4, 'off')
            end
            out.v1Acorr = yProfile1;
            out.y1Acorr = yshift1;
            out.acorrFit1 = gaussFit1;

            if obj.display
                val = feval(gaussFit2, -p.maxshift:p.maxshift);
    %             i1=find(yProfile>yMax_val/2,1,'first');
    %             i2=find(yProfile>yMax_val/2,1,'last');
    %             fwhm = (i2-i1)*maxprecision;
    %             size_oneK = fwhm - 2.355*mean(locs1.locprecnm);
                ax5=obj.setoutput('fwhm_auto_l2');
                plot(ax5, yshift2, yProfile2,'r');
                hold(ax5, 'on')
    %             yline(ax4, yMax_val/2)
                plot(ax5, -p.maxshift:p.maxshift, val);
                yline(ax5, bg2);
    %             title(ax4, ['FWHM = ' num2str(fwhm) ' nm; width_k = ' num2str(size_oneK)])
                xlabel(ax5, 'y shift (nm)')
                ylabel(ax5, 'auto corr')
                hold(ax5, 'off')
            end
            out.v2Acorr = yProfile2;
            out.y2Acorr = yshift2;
            out.acorrFit2 = gaussFit2;
%             out.fwhm1 = fwhm;
%             out.size_oneK = size_oneK;

        end
        function pard=guidef(obj)
            pard=guidef;
        end

    end
end

function [yshift, yProfile, gaussFit, bg] = auto_corr_analysis(sxf, p, maxprecision)
    sxhr=imresize(sxf,p.pixrec/maxprecision,'cubic');
    [~,linind]=max(sxhr(:));
    [xm,ym]=ind2sub(size(sxhr),linind);
% out.dxline=(xm-size(sx12hr,1)/2)*maxprecision;
    yProfile = sxhr(xm,:);
    yProfile = yProfile-min(yProfile); 
%             yMax_val = max(yProfile);
%     gauss2 = 'a1*exp(-((x-b)/c1)^2) + a2*exp(-((x-b)/c2)^2)+d';
    gauss2 = 'a1*exp(-((x-b)/c1)^2) + a2*exp(-((x-b)/c2)^2)';
    yshift = (-size(sxhr,1)/2:size(sxhr,1)/2)*maxprecision;
    yshift = yshift(1:end-1);
%     gaussFit = fit(yshift',yProfile', gauss2, 'StartPoint',[max(yProfile)*0.9 max(yProfile)*0.1 0 10 100 0]);
    gaussFit = fit(yshift',yProfile', gauss2, 'StartPoint',[max(yProfile)*0.9 max(yProfile)*0.1 0 25 100]);
%     bg = gaussFit.a2+gaussFit.d;
    bg = gaussFit.a2;
end
function pard=guidef
pard.inputParameters={'se_sitepixelsize','numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi','layer2_'};

pard.lineselect.object=struct('Style','popupmenu','String',{{'line1','line2','rotationline'}});
pard.lineselect.position=[1,1];
pard.lineselect.Width=4;

%maximum shift
% fit window

pard.maxshiftt.object=struct('Style','text','String','max shift nm');
pard.maxshiftt.position=[2,1];
pard.maxshiftt.Width=1;

pard.maxshift.object=struct('Style','edit','String','50');
pard.maxshift.position=[2,2];
pard.maxshift.Width=.5;

pard.fitwint.object=struct('Style','text','String','filter nm');
pard.fitwint.position=[2,3];
pard.fitwint.Width=1;

pard.fitsig.object=struct('Style','edit','String','5');
pard.fitsig.position=[2,4];
pard.fitsig.Width=.5;

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


