classdef lineprofile<interfaces.SEEvaluationProcessor
%     Calculates the line-profile along user-defined direction and fits it
%     with a Gaussian or double-Gaussian model
    properties
        boundary
    end
    methods
        function obj=lineprofile(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj, p)
            % Import the kimograph
%             site = obj(1).locData
%             runNPC3Dfitting(obj,inp)
            layers=find(p.sr_layerson);
            if obj.display
                ax=obj.setoutput('profile');
                ax.Position(3)=0.6;
                hold(ax,'off')
                ax4=obj.setoutput('zprofile');
                hold(ax4,'off')
            end
            t1={};
            lt={};
            if p.setbinwidth
                pixelsize=p.binwidth;
            else
                pixelsize=p.se_sitepixelsize;
            end
            for k=1:length(layers)
                locs=obj.getLocs({'xnmrot','ynmrot','znm','locprecnm','locprecznm'},'layer',layers(k),'size',p.se_siteroi/2);  
                nbins=-200:pixelsize:200;
                hc=histcounts(locs.ynmrot,nbins);
                nbinsf=nbins(1:end-1)+(nbins(2)-nbins(1))/2;
%                 mp=median(locs.ynmrot);
%                 sd=std(locs.ynmrot);
% %                 d=30;
%                 amp=max(hc);
%                 ft=fittype( @(a1,b1,s,a2,b2,x) a1*exp(-(x-b1).^2/2/s^2)+a2*exp(-(x-b2).^2/2/s^2));
% %                 ft=fittype( @(a1,b1,s,a2,b2,c,x) a1*normcdf((x-b1)/s)+a2*normcdf((x-b2)/s)+c);
% %                 hc=cumsum(hc);
%                 dstart=50;
%                 fitp=fit(nbins(1:end-1)',hc',ft,'StartPoint',[amp mp-dstart/2 sd amp mp+dstart/2],'Lower',[0 -inf 0 0 -inf]);
                [fitp,fitprof,fittext]=fitgeneralprofile(hc,nbinsf,p,mean(locs.locprecnm));
                if obj.display
                    plot(ax,nbins(1:end-1),hc,nbinsf,fitprof);
                    hold(ax,'on');
                    t1=[t1 {['layer ' num2str(k) ':']} fittext];
                    lt=[lt {['layer ' num2str(k) ':']} 'fit'];
    %                 dist=abs(fitp.b2-fitp.b1);
    %                 title(ax,dist)
                end
                if ~isempty(locs.znm)
                    % 201412: Yu-Le added this part for z-profile
                    % copied from make_lineprofiles
                    binwidth = pixelsize;
                    
                    z=double(locs.znm);
                    locprecznm=double(locs.locprecznm);
                    if isempty(z)
                        z=x*0;
                        locprecznm=z;
                    end

                    minzh=max(-750,min(z));
                    maxzh=min(750,max(z));
                    n=minzh-3*binwidth:binwidth:maxzh+3*binwidth;
                    profz=hist(z,n);profz([1 end])=[];n([1 end])=[];
                    
                    fwhm=getFWHM(profz,n);
%                     t3{end+1}=['FWHM: ' 9  num2str(fwhm)];
                    [fitp,fitprof,fittxt]=fitgeneralprofile(profz,n,p,fwhm/2.6);
                    if obj.display
                        axes(ax4)
                        plot(ax4,n,profz);
                        hold(ax4,'on')
                        xlabel(ax4,'z (nm)')
                        ylabel(ax4,'counts')
                        plot(ax4,n,fitprof,'k--')
                    end
%                     t3(end+1:end+length(fittxt))=fittxt;

%                     axes(ax5)
%                     plot(x,z,'.')
%                     xlabel(ax5,'Position along line ROI (nm)')
%                     ylabel(ax5,'z (nm)')
%                     axis equal tight
%                     hold on
                end
                out.fitp{k}=fitp;
            end
            if obj.display
                legend(ax,lt)
                pos=[.75,0.025,.25,.95];
                fontsize=12;
                uicontrol('Parent',ax.Parent,'style','text','String',t1,'Units','normalized','Position',pos,'FontSize',fontsize,'HorizontalAlignment','left')  
                out.model=p.fitmodel.selection;
            end
%             out.distance=dist;
%             out.std=fitp.s;
        end
     
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end

end



function pard=guidef(obj)
pard.fitmodel.object=struct('Style','popupmenu','String',{{'Gaussian','Step','Disk','Ring','Double Gaussian'}});
pard.fitmodel.position=[1,1];
pard.fitmodel.Width=2;

pard.restrictsigma.object=struct('String','sigma=<locp>','Style','checkbox');
pard.restrictsigma.position=[3,1];
pard.restrictsigma.Width=4;
 p(1).value=0; p(1).on={}; p(1).off={'binwidth'};
p(2).value=1; p(2).on={'binwidth'}; p(2).off={};
pard.setbinwidth.object=struct('String','set binwidth (nm):','Style','checkbox','Value',1,'Callback',{{@obj.switchvisible,p}});
pard.setbinwidth.position=[2,1];
pard.setbinwidth.Width=3;
pard.setbinwidth.Tooltip='Set bin width for profile. If not set, use pixel size from ROI manager';
pard.binwidth.object=struct('String','2','Style','edit');
pard.binwidth.position=[2,3.5];
pard.binwidth.Width=0.5;
pard.binwidth.TooltipString='Binwidth for profiles. If not checked, use pixel size of reconstruction';
pard.setbinwidth.TooltipString=pard.binwidth.TooltipString;

% pard.dxt.Width=3;
pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi','se_sitepixelsize'};
pard.plugininfo.type='ROI_Evaluate';
pard.plugininfo.description='Calculates the line-profile along user-defined direction and fits it with a Gaussian or double-Gaussian model';
end

function M = calMeasurement(x,y,qx,qy,Size)
    ref = interp1(x,y,qx);
    rightIdx = qy > ref;
    leftIdx = ~rightIdx;
    leftA = sum(y);
    rightA = Size(1)*Size(2)-sum(y);
    leftD = sum(leftIdx)/leftA;
    rightD = sum(rightIdx)/rightA;
    M = 0-((leftD-1)^2 + (rightD-0)^2)^(1/2);
end

function [fwhm,fwhmind]=getFWHM(profile,x)
    [mp, ip]=max(profile);
    i1=find(profile(1:ip)>mp/2,1,'first');
    i2=find(profile(ip:end)>mp/2,1,'last')+ip-1;
    if isempty(i2)||isempty(i1)
        fwhm=[];
        fwhmind=1;
    else
        if i1==i2
            i1=i1-1;i2=i2+1;
        end
        
    fwhm=x(i2)-x(i1);
    fwhmind=i2-i1;
    end
end
