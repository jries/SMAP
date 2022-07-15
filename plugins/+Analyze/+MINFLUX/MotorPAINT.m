classdef MotorPAINT<interfaces.DialogProcessor
    % LINEPROFILE Calculates profiles along a linear ROI and fits it with a
    % model of choice. Flat: step function convolved with Gaussian
    % (=Erf). Disk: Projection of a homogeneously filled disk, convolved
    % with Gaussian. Ring: Projection of a ring, convolved with
    % Gaussian. Distance: Two Gaussians in a distance d.
    methods
        function obj=MotorPAINT(varargin)        
            obj@interfaces.DialogProcessor(varargin{:});
            obj.showresults=true;
        end
        
        function out=run(obj,p)  
            out=[];

            [locs,indin]=obj.locData.getloc({'numberInGroup','groupindex','xnm','ynm','znm','time','frame'},'layer',1,'Position','fov');
            findin=find(indin);
            gn=unique(locs.groupindex(locs.numberInGroup>=p.minlen));
            isz=~isempty(locs.znm);
            
            locsout.xnm=[];
            locsout.ynm=[];
            
            locsout.frame=[];
            locsout.tracklength=[];locsout.tracknumber=[];locsout.trackangle=[];
            if isz
                locsout.znm=[];
            end
            tind=1;
            t=tic;
            
            axall=obj.initaxis('raw');
            hold(axall,'off');
            if ~isfield(obj.locData.loc,'trackangle')
                obj.locData.setloc('trackangle',0*obj.locData.loc.xnm)
                obj.locData.setloc('tracklength',0*obj.locData.loc.xnm)
            end
            for k=1:length(gn)
                indh=locs.groupindex==gn(k);

                if p.skipfirst>0
                    outn=find(indh,p.skipfirst,'first');
                    indh(outn)=false;
                end
                if sum(indh)<p.minlen
                    continue
                end
                
                x=locs.xnm(indh);
                y=locs.ynm(indh);
                frame1=locs.frame(find(indh,1,'first'));


                lennm=sqrt(diff(quantile(x,[.1,0.9]))^2+diff(quantile(y,[.1,0.9]))^2);

%                 lennm=sqrt((x(end)-x(1))^2+(y(end)-y(1))^2);
                if lennm<p.minlennm
                    continue
                end
                plot(axall,x,y);
                hold(axall,'on');
                time=locs.time(indh);
                
                %p. splitmerge (T), splitmergestep([]), stepfunction (mean), coarsneness
                p.splitmerge=~isempty(p.splitmergestep);
                p.stepfunction=p.stepfunctionm.selection; 

                [xr,yr,angle]=rotateCenterCoordinates(x,y,time);
                switch p.smoothingfunction.selection
                    case 'simple smooth'
                        istep=smoothtrackind(xr,p);
                    case 'step finder'
                         istep=findstepsMINFLUX(xr,p);
                end
%                 try
               
                
%                 catch err
%                     err
%                     rethrow(err)
%                     continue;
%                 end
                [xs,locsstep]=stepvalue(x,istep);
                
                ys=stepvalue(y,istep);
                
                tracklength=(sqrt((xs(end)-xs(1))^2+(ys(end)-ys(1))^2));
                
                locsout.xnm(end+1:end+length(xs))=single(xs);
                locsout.ynm(end+1:end+length(xs))=single(ys);
                locsout.frame(end+1:end+length(xs))=single(time(istep));
                locsout.tracklength(end+1:end+length(xs))=single(ones(size(xs))*tracklength);
                locsout.tracknumber(end+1:end+length(xs))=single(ones(size(xs))*gn(k));
                locsout.trackangle(end+1:end+length(xs))=single(ones(size(xs))*angle);
                if isz
                    z=locs.znm(indh);
                    zs=stepvalue(z,istep);
                    locsout.znm(end+1:end+length(zs))=single(zs);
                    zsall{tind}=zs;
                    zrawall{tind}=z;
                end

                xsall{tind}=xs;ysall{tind}=ys;
                xrawall{tind}=x;yrawall{tind}=y;
                angleall(tind)=angle;
                frameall(tind)=frame1;
                tind=tind+1;

                obj.locData.loc.trackangle(findin(indh))=angle;
                obj.locData.loc.tracklength(findin(indh))=tracklength;

                if toc(t)>1
                    obj.status(['MINFLUX tracks ' num2str(k/length(gn)*100,'%2.1f') '%, ' num2str(k) ' / ' num2str(length(gn))]); drawnow
                    t=tic;
                end
            end
            axis(axall,"ij")            
            axis(axall,"equal")

            [~, filename]=fileparts(obj.locData.files.file(1).name);
            obj.locData.addfile(['tracks_s' num2str(p.splitmergestep) '_' num2str(obj.locData.files.filenumberEnd) '_' filename]);
            % obj.locData.files.file(end).info.simulationParameters=obj.getGuiParameters;
            if p.addtracks
            obj.locData.addLocData(locsout);

            obj.locData.filter
            %            try
            initGuiAfterLoad(obj);
            obj.locData.filter
            end

            ax=obj.initaxis('tracks');
            hold(ax,'off')
            for k=1:length(xsall)
                plot(ax,xsall{k},ysall{k})
                hold(ax,'on')
            end
            axis(ax,"equal")
            axis(ax,"ij")
            ax.XLim=axall.XLim;
            ax.YLim=axall.YLim;

            axc=obj.initaxis('angle');
            cmap=hsv(256);
            hold(axc,'off')
            for k=1:length(xsall)
                cind=max(1,ceil(angleall(k)/2/pi*255));
                plot(axc,xsall{k},ysall{k},'Color',cmap(cind,:))
                hold(axc,'on')
            end
            axis(axc,"equal")
            axis(axc,"ij")
            axc.XLim=axall.XLim;
            axc.YLim=axall.YLim;
            if isz
                r=rand(length(xsall),1);
                [~,p.indsort]=sort(r);

                ax3Dt=obj.initaxis('track 3D');
                p.plotsides=true;
                p.use=[];
                plot3Dtracks(ax3Dt,xsall,ysall,zsall,[-100 100 -100],p)
                
                p.use=frameall>157000 & frameall<213000;
                p.plotsides=false;
                ax3Dtr=obj.initaxis('track 3D raw');
                plot3Dtracks(ax3Dtr,xrawall,yrawall,zrawall,[-100 100 -100],p)
%                 hold(ax3Dtr,'off')
%                 skipfirst=10;
%                 c=jet(length(xrawall));
%                 for tr=1:length(xrawall)
%                     xh=xrawall{tr}(skipfirst+1:end);
%                     yh=yrawall{tr}(skipfirst+1:end);
%                     zh=zrawall{tr}(skipfirst+1:end);
% 
%                     hl=plot3(ax3Dtr,xh,yh,zh,'Color',c(tr,:),'LineWidth',.5);
%                     hold(ax3Dtr,'on')
% 
%                 end
%                 axis(ax3Dtr,'tight')
%                 axis(ax3Dtr,'equal');
%                 view(ax3Dtr,70,27)
                

            end

        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end




function pard=guidef(obj)

%  p(1).value=0; p(1).on={}; p(1).off={'binwidth'};
% p(2).value=1; p(2).on={'binwidth'}; p(2).off={};
% pard.setbinwidth.object=struct('String','set binwidth (nm) (otherwise: pixelsize):','Style','checkbox','Value',1,'Callback',{{@obj.switchvisible,p}});
% pard.setbinwidth.position=[1,1];
% pard.setbinwidth.Width=3;
% 
pard.minlent.object=struct('String','minimum track length (locs)','Style','text');
pard.minlent.position=[1,1];
pard.minlent.Width=1.5;
pard.minlen.object=struct('String','100','Style','edit');
pard.minlen.position=[1,2.5];
pard.minlen.Width=0.5;

pard.minlennmt.object=struct('String','minimum track length (nm)','Style','text');
pard.minlennmt.position=[1,3];
pard.minlennmt.Width=1.5;
pard.minlennm.object=struct('String','100','Style','edit');
pard.minlennm.position=[1,4.5];
pard.minlennm.Width=0.5;

pard.splitmergestept.object=struct('String','step size (nm)','Style','text');
pard.splitmergestept.position=[2,1];
pard.splitmergestept.Width=1.5;
pard.splitmergestep.object=struct('String','8','Style','edit');
pard.splitmergestep.position=[2,2.5];
pard.splitmergestep.Width=0.5;


pard.smoothingfunction.object=struct('String',{{'step finder','simple smooth'}},'Style','popupmenu');
pard.smoothingfunction.position=[2,3];
pard.smoothingfunction.Width=2;

pard.coarsenesst.object=struct('String','Coarseness','Style','text');
pard.coarsenesst.position=[3,1];
pard.coarsenesst.Width=1.5;
pard.coarseness.object=struct('String','1','Style','edit');
pard.coarseness.position=[3,2.5];
pard.coarseness.Width=0.5;

pard.stepfunctiont.object=struct('String','Function','Style','text');
pard.stepfunctiont.position=[3,3];
pard.stepfunctiont.Width=1.;
pard.stepfunctionm.object=struct('String',{{'mean','median'}},'Style','popupmenu');
pard.stepfunctionm.position=[3,4.];
pard.stepfunctionm.Width=1;

pard.skipfirstt.object=struct('String','Skip first:','Style','text');
pard.skipfirstt.position=[4,1];
pard.skipfirstt.Width=1.5;
pard.skipfirst.object=struct('String','0','Style','edit');
pard.skipfirst.position=[4,2.5];
pard.skipfirst.Width=0.5;

pard.addtracks.object=struct('String','Add tracks','Style','checkbox','Value',1);
pard.addtracks.position=[5,1];
pard.addtracks.Width=2;

pard.plugininfo.description=sprintf('');
pard.plugininfo.type='ProcessorPlugin';
end