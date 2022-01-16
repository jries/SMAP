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

            locs=obj.locData.getloc({'numberInGroup','groupindex','xnm','ynm','time'},'layer',1,'Position','fov');
            
            gn=unique(locs.groupindex(locs.numberInGroup>=p.minlen));
            
            locsout.xnm=[];
            locsout.ynm=[];
            locsout.frame=[];
            locsout.tracklength=[];locsout.tracknumber=[];
            tind=1;
            t=tic;
            
            for k=1:length(gn)
                indh=locs.groupindex==gn(k);
                
                x=locs.xnm(indh);
                y=locs.ynm(indh);

                lennm=sqrt(diff(quantile(x,[.1,0.9]))^2+diff(quantile(y,[.1,0.9]))^2);

%                 lennm=sqrt((x(end)-x(1))^2+(y(end)-y(1))^2);
                if lennm<p.minlennm
                    continue
                end
                time=locs.time(indh);
                
                %p. splitmerge (T), splitmergestep([]), stepfunction (mean), coarsneness
                p.splitmerge=~isempty(p.splitmergestep);
                p.stepfunction=p.stepfunctionm.selection; 

                [xr,yr]=rotateCenterCoordinates(x,y,time);
%                 try
                istep=findstepsMINFLUX(xr,p);
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
                xsall{tind}=xs;ysall{tind}=ys;tind=tind+1;

                if toc(t)>1
                    obj.status(['MINFLUX tracks ' num2str(k/length(gn)*100,'%2.1f') '%, ' num2str(k) ' / ' num2str(length(gn))]); drawnow
                    t=tic;
                end
            end
            
            [~, filename]=fileparts(obj.locData.files.file(1).name);
            obj.locData.addfile(['tracks_s' num2str(p.splitmergestep) '_' num2str(obj.locData.files.filenumberEnd) '_' filename]);
            % obj.locData.files.file(end).info.simulationParameters=obj.getGuiParameters;
            obj.locData.addLocData(locsout);

            obj.locData.filter
            %            try
            initGuiAfterLoad(obj);
            obj.locData.filter

            ax=obj.initaxis('tracks');
            hold(ax,'off')
            for k=1:length(gn)
                plot(ax,xsall{k},ysall{k})
                hold(ax,'on')
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


pard.plugininfo.description=sprintf('');
pard.plugininfo.type='ProcessorPlugin';
end