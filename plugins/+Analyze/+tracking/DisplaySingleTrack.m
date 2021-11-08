classdef DisplaySingleTrack<interfaces.DialogProcessor
    %Interactive analysis of SPT data
    methods
        function obj=DisplaySingleTrack(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.inputParameters={'sr_pixrec','numberOfLayers','sr_pos','sr_size','layers','sr_layerson','cam_pixelsize_nm'};
            obj.showresults=true;
        end
        function out=run(obj,p)
            out=[];
            locs=obj.locData.getloc({'xnm','ynm','znm','frame','time','xnmline','ynmline'},'layer',find(obj.getPar('sr_layerson')), 'grouping','ungrouped','Position','roi');
            if ~isempty(locs.time)
                time=locs.time;
                xax='time (ms)';
                ds=quantile(diff(locs.time),0.05);
                
            else
                time=locs.frame;
                xax='frame (xx)';
                ds=quantile(diff(locs.frame),0.05);
            end
            if ~isempty(locs.xnmline)
                x=locs.xnmline;
                y=locs.ynmline;
            else
                x=locs.xnm;
                y=locs.ynm;
            end
            
            mintimestep=100; 
            %keep longest continuous track
            indgood=true(size(time));
            dtmax=inf;
            while dtmax>mintimestep
                dt=diff(time(indgood));
                findgood=find(indgood);
                [dtmax,indm]=max(dt);
                if dtmax>mintimestep
                    if indm>length(dt)/2
                        indgood(findgood(indm):end)=false;
                    else
                        indgood(1:findgood(indm))=false;
                    end
                end
            end
            
            x=x(indgood);
            time=time(indgood);
            y=y(indgood);
            
            median(diff(time))
            
            ax=obj.initaxis('correlation');
            h=histcounts(x,min(x):1:max(x));
            xc=myxcorr(h,h);
            plot(ax,xc);
            xlabel('delta x (nm)');
            ylabel('auto corr')
            
            if 0% p.stepsize>0
                ax=obj.initaxis('steppos');
                xm=mod(x,p.stepsize);
                hxm=histcounts(xm,0:p.stepsize);
                plot(hxm);
                [~,stepshift]=max(hxm);



                axx=obj.initaxis('x');
                steps=(ceil(min(x)/p.stepsize)*p.stepsize:p.stepsize:max(x)+stepshift)-stepshift;
                stepsa=vertcat(steps,steps);
                tsteps=vertcat(ones(size(steps))*min(time),ones(size(steps))*max(time));
                hold(axx,'off')
                plot(axx,tsteps,stepsa,'c')
                hold(axx,'on')
            
            else
                axx=obj.initaxis('x');
                hold(axx,'off')
            end
                
            timenorm=true;
            if timenorm
                timeplot=(time-time(1))/1000;
                xax='time (s)';
            else
                timeplot=time;
            end
            plot(axx,timeplot,x,'m');
            hold(axx,'on')
            xlabel(axx,xax)
            ylabel(axx,'position (nm)')
            windowsize=p.filtersize*ds;
            xf=runningWindowAnalysis(time,x,time,windowsize,p.filtermode.selection);
            
            plot(axx,timeplot,xf,'b');
            
%             plot(ax,time,y,'m');
            
            
            try
            %step finder
            ax=obj.initaxis('stepfind');
            stepsize=p.stepsize;

            inds=findchangepts(x,'MaxNumChanges',round((max(x)-min(x))/max(20,stepsize)));

            inds=[0 ;inds ;length(x)+1];
            mv=zeros(length(inds)-1,1);
            tv=mv;
            for k=1:length(inds)-1
                mv(k)=mean(x(inds(k)+1:inds(k+1)-1));
                tv(k)=timeplot(inds(k)+1);

            end
            stairs(axx,tv,mv,'k')
            
            
            dmv=diff(mv);
            for k=1:length(dmv)
                text(axx,tv(k+1),mean(mv(k:k+1)),num2str(dmv(k),'%2.0f'))
            end
            histogram(ax,dmv,min(dmv)-5:5:max(dmv)+5)
            
            ax=obj.initaxis('x-y');
            plot(ax,x,y,'.-')
            
            obj.initaxis('x','keep');
            catch err
                disp('Signal processing Toolbox needed')
            end
            if p.makemovie
                ts=min(time):p.frametime:max(time);
                f=figure(99);
                f.Position(1)=1;f.Position(3)=1280;
                ax=gca;

                delete(ax.Children)
                
                axis(ax,'equal');
%                 axis(ax,'off');
                    xlim(ax,[min(x)-10 max(x)+10])
                    ylim(ax,[min(y)-10 max(y)+10])
                    hold(ax,'on')
                    plot([min(x)-5 min(x)+10-5], [min(y)-5 min(y)-5],'k','LineWidth',3)
                    ax.XTick=[];
                    ax.YTick=[];
                for k=1:length(ts)
                    indh=time<=ts(k);
                    xh=x(indh);
                    yh=y(indh);
                    th=time(indh);
%                     hold(ax,'off')
                    hl=plot(ax,xh,yh,'k');
                    hold(ax,'on')
                    hd=plot(ax,xh(end),yh(end),'ro','MarkerFaceColor','r','MarkerSize',10);
                    hb=plot(ax,xh,yh,'b.');
                    tpassed=ts(k)-ts(1);
                    ht=text(ax,double(min(x)),double(max(y)),[num2str(tpassed,'%3.0f') ' ms'],'FontSize',15);
                    
                    drawnow
                    Fr(k)=getframe(ax);
                    delete(hd)
                    delete(hl)
                    delete(ht)
                    delete(hb)
                end
                smlfile=obj.getPar('lastSMLFile');
                if ~isempty(smlfile)
                    pfad=fileparts(smlfile);
                else
                    pfad=fileparts(obj.locData.files.file(1).name);
                end
%                 fo=strrep(fo,'_sml.mat','.mp4');
                [file,pfad]=uiputfile([pfad filesep '*.mp4']);
                if file
                    mysavemovie(Fr,[pfad  file],'FrameRate',20)
                end
                  
                    
                   
            end

        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end


function pard=guidef

pard.filtert.object=struct('String','filter size','Style','text');
pard.filtert.position=[1,1];

pard.filtersize.object=struct('String','10','Style','edit');
pard.filtersize.position=[1,2];
pard.filtermode.object=struct('String',{{'median','mean'}},'Style','popupmenu');
pard.filtermode.position=[1,3];


pard.stepsizet.object=struct('String','Step size (nm)','Style','text');
pard.stepsizet.position=[2,1];
pard.stepsize.object=struct('String','0','Style','edit');
pard.stepsize.position=[2,2];


pard.makemovie.object=struct('String','Make movie','Style','checkbox');
pard.makemovie.position=[3,1];

pard.frametimet.object=struct('String','frame time (ms)','Style','text');
pard.frametimet.position=[3,2];

pard.frametime.object=struct('String','10','Style','edit');
pard.frametime.position=[3,3];
% pard.analysismode.Width=3;
% % ,'grid based diffusion coefficients'
% 
% pard.lentrackst.object=struct('String','Minimum length of tracks','Style','text');
% pard.lentrackst.position=[2,1];
% pard.lentrackst.Width=2.5;
% pard.lentrackst.TooltipString=sprintf(['set this keyword to eliminate all trajectories with \n'...
%             ' fewer than param.good valid positions.  This is useful \n'...
%            'due to blinking noise particles in the data stream.']);
% pard.lentracks.object=struct('String','2','Style','edit');
% pard.lentracks.position=[2,3.5];
% pard.lentracks.Width=.5;
% pard.lentracks.TooltipString=pard.lentrackst.TooltipString;
% 
% pard.cutoffimmobilet.object=struct('String','Cutoff immobile log10(D(um^2/s))','Style','text');
% pard.cutoffimmobilet.position=[3,1];
% pard.cutoffimmobilet.Width=2.5;
% pard.cutoffimmobilet.TooltipString=sprintf(['set this keyword to eliminate all trajectories with \n'...
%             ' fewer than param.good valid positions.  This is useful \n'...
%            'due to blinking noise particles in the data stream.']);
% pard.cutoffimmobile.object=struct('String','-3','Style','edit');
% pard.cutoffimmobile.position=[3,3.5];
% pard.cutoffimmobile.Width=.5;
% pard.cutoffimmobile.TooltipString=pard.cutoffimmobilet.TooltipString;
% 
% pard.saveD.object=struct('String','Write diffusion coefficient to localization data','Style','checkbox');
% pard.saveD.position=[4,1];
% pard.saveD.Width=3;


pard.plugininfo.description=sprintf('Interactive analysis of SPT data');
pard.plugininfo.type='ProcessorPlugin';
% pard.plugininfo.Name='AnalyzeSPT';
end

