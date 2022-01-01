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

            if p.usefilter
                windowsize=p.filtersize*ds;
                xf=runningWindowAnalysis(time,x,time,windowsize,p.filtermode.selection);            
                plot(axx,timeplot,xf,'b');
            end
            
%             plot(ax,time,y,'m');
            
            
%             try
            %step finder
            ax=obj.initaxis('stepfind');
%             stepsize=p.stepsize;
            p.fitmean=contains(p.fitmode.selection,'mean');
            steps=AutoStepfinderRies(x,p);

%             inds=findchangepts(x,'MaxNumChanges',round((max(x)-min(x))/max(20,stepsize)));

            inds=[1 ;steps.properties.IndexStep ];
            mv=[steps.properties.LevelBefore(1) ;steps.properties.LevelAfter];
            tv=timeplot(inds);
%             mv=zeros(length(inds)-1,1);
%             tv=mv;
%             for k=1:length(inds)-1
%                 mv(k)=mean(x(inds(k)+1:inds(k+1)-1));
%                 tv(k)=timeplot(inds(k)+1);
% 
%             end
            stairs(axx,tv,mv,'k')
       
            
            dmv=diff(mv);
            for k=1:length(dmv)
                text(axx,tv(k+1),mean(mv(k:k+1)),num2str(dmv(k),'%2.0f'))
            end
            histogram(ax,dmv,min(dmv)-5:5:max(dmv)+5)
%             catch err
%                 disp('Signal processing Toolbox needed')
%             end              
            ax=obj.initaxis('x-y');
            plot(ax,x,y,'.-')
            axis(ax,'equal')
            
            obj.initaxis('x','keep');

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
            pard=guidef(obj);
        end
    end
end


function pard=guidef(obj)

p(1).value=0; p(1).on={}; p(1).off={'filtert','filtersize','filtermode'};
p(2).value=1; p(2).on=p(1).off; p(2).off={};
pard.usefilter.object=struct('String','filter','Style','checkbox','Callback',{{@obj.switchvisible,p}},'Value',1);
pard.usefilter.position=[1,1];
pard.usefilter.Width=0.5;
pard.filtert.object=struct('String','kernel','Style','text');
pard.filtert.position=[1,1.5];
pard.filtert.Width=0.5;
pard.filtersize.object=struct('String','10','Style','edit');
pard.filtersize.position=[1,2];
pard.filtersize.Width=0.5;
pard.filtermode.object=struct('String',{{'median','mean'}},'Style','popupmenu');
pard.filtermode.position=[1,2.5];


% pard.stepsizet.object=struct('String','Step size (nm)','Style','text');
% pard.stepsizet.position=[2,1];
% pard.stepsize.object=struct('String','0','Style','edit');
% pard.stepsize.position=[2,2];

p(1).value=0; p(1).on={}; p(1).off={'frametimet','frametime'};
p(2).value=1; p(2).on=p(1).off; p(2).off={};
pard.makemovie.object=struct('String','Make movie','Style','checkbox','Callback',{{@obj.switchvisible,p}});
pard.makemovie.position=[2,1];

pard.frametimet.object=struct('String','frame time (ms)','Style','text','Visible','off');
pard.frametimet.position=[2,2];
pard.frametime.object=struct('String','10','Style','edit','Visible','off');
pard.frametime.position=[2,3];

pard.overshoott.object=struct('String','Coarsness','Style','text');
pard.overshoott.position=[3,1];
pard.overshoot.object=struct('String','1','Style','edit');
pard.overshoot.position=[3,2];
pard.overshoot.Width=0.5;


pard.fitmodet.object=struct('String','fit using','Style','text');
pard.fitmodet.position=[3,3];
pard.fitmode.object=struct('String',{{'mean','median'}},'Style','popupmenu');
pard.fitmode.position=[3,4];

% pard.SMaxTresholdt.object=struct('String','Fitting threshold','Style','text');
% pard.SMaxTresholdt.position=[3,3];
% pard.SMaxTreshold.object=struct('String','0.1','Style','edit');
% pard.SMaxTreshold.position=[3,4];
% pard.SMaxTreshold.Width=0.5;

p(1).value=0; p(1).on={}; p(1).off={'manualmodestepst','manualmodesteps'};
p(2).value=1; p(2).on=p(1).off; p(2).off={};
pard.manualon.object=struct('String','Manual mode','Style','checkbox','Callback',{{@obj.switchvisible,p}});
pard.manualon.position=[4,1];
pard.manualmodestepst.object=struct('String','Manual steps','Style','text','Visible','off');
pard.manualmodestepst.position=[4,2];
pard.manualmodesteps.object=struct('String','10','Style','edit','Visible','off');
pard.manualmodesteps.position=[4,3];
pard.manualmodesteps.Width=0.5;

% p(1).value=0; p(1).on={}; p(1).off={'meanbaset','meanbase'};
% p(2).value=1; p(2).on=p(1).off; p(2).off={};
% pard.PostProcessOn.object=struct('String','Post processing','Style','checkbox','Callback',{{@obj.switchvisible,p}});
% pard.PostProcessOn.position=[5,1];
% pard.meanbaset.object=struct('String','Threshold','Style','text','Visible','off');
% pard.meanbaset.position=[5,2];
% pard.meanbaset.Width=0.5;
% pard.meanbase.object=struct('String','0','Style','edit','Visible','off');
% pard.meanbase.position=[5,2.5];
% pard.meanbase.Width=0.5;

p(1).value=0; p(1).on={}; p(1).off={'noisemaxdistt','noisemaxdist'};
p(2).value=1; p(2).on=p(1).off; p(2).off={};
pard.noseeston.object=struct('String','Noise estimation','Style','checkbox','Callback',{{@obj.switchvisible,p}});
pard.noseeston.position=[5,3];
pard.noseeston.Width=1;
pard.noisemaxdistt.object=struct('String','Range','Style','text','Visible','off');
pard.noisemaxdistt.position=[5,4];
pard.noisemaxdistt.Width=0.5;
pard.noisemaxdist.object=struct('String','100','Style','edit','Visible','off');
pard.noisemaxdist.position=[5,4.5];
pard.noisemaxdist.Width=0.5;
%post processing threshold  ??
% noise estimation range max_range?

pard.plugininfo.description=sprintf('Interactive analysis of SPT data');
pard.plugininfo.type='ProcessorPlugin';
% pard.plugininfo.Name='AnalyzeSPT';
end

