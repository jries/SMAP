classdef StatsVsTime<interfaces.DialogProcessor
    % StatsVsTime calculates localization statistics in dependence of the
    % filenumber or of the frame
    methods
        function obj=StatsVsTime(varargin)           
            obj@interfaces.DialogProcessor(varargin{:}) 
            obj.inputParameters={'sr_layerson','layernames'};
            obj.showresults=true;
        end
        
        function out=run(obj,p)
            out=[];
           
            fields={'filenumber','frame','phot','locprecnm','PSFxnm','numberInGroup','bg'};
            if ~isfield(obj.locData.loc,'bg')
                obj.locData.setloc('bg',obj.locData.loc.bg1)
                disp('bg1 used for bg')
                obj.locData.regroup;
            end
            if isfield(obj.locData.loc,'znm')
                fields(end+1:end+2)={'znm','locprecznm'};
            end
            if p.useroi
                position='roi';
            else
                position='all';
            end
            usefields={{'photons','Nloc'},{'photons','mu'},{'lifetime','mu'},{'background','max'},{'PSFxnm','max'},{'frames','falloff'},{'locprec','max'},{'locprec','rising'}};
            layers=find(p.sr_layerson);            
            
            switch p.timefield.selection
                case 'files'
                    maxfn=max(obj.locData.getloc('filenumber').filenumber);
                    for filen=1:maxfn
                        if p.filter
                            for m=length(layers):-1:1
                                locs{m}=obj.locData.getloc(fields,'layer',layers(m),'position',position,'filenumber',filen);
                                modetxt{m}=['layer' num2str(layers(m))];
                            end
                        else
                            locs{2}=obj.locData.getloc(fields,'position',position,'grouping','grouped','filenumber',filen);
                            locs{1}=obj.locData.getloc(fields,'position',position,'grouping','ungrouped','filenumber',filen);
                            modetxt{2}='grouped';
                            modetxt{1}='ungrouped';
                        end
                        
                        stats=make_statistics2(locs,p,false);
                        for f=1:length(usefields)
                            fh=usefields{f};
                            if isfield(stats,fh{1})
                            stath=stats.(fh{1}).(fh{2});
                            for l=1:length(stath)
                                ps.(fh{1}).(fh{2}){l}(filen)=stath(l);
                            end
                            end
                        end
                        
                    end
                    filns=1:maxfn;
                   xl='filenumber';

                case 'frames'                   
                    if p.filter
                        for m=length(layers):-1:1
                            locs{m}=obj.locData.getloc(fields,'layer',layers(m),'position',position);
                            modetxt{m}=['layer' num2str(layers(m))];
                        end
                    else
                        locs{2}=obj.locData.getloc(fields,'position',position,'grouping','grouped');
                        locs{1}=obj.locData.getloc(fields,'position',position,'grouping','ungrouped');
                        modetxt{2}='grouped';
                        modetxt{1}='ungrouped';
                    end
                    maxf=0;
                    for k=1:length(locs)
                        maxf=max(maxf,max(locs{k}.frame));
                    end
                    df=ceil(maxf/p.framewindows);
                    frames=[1:df:maxf maxf];
                    for f=1:length(frames)-1
                        for k=1:length(locs)
                            indf=locs{k}.frame>=frames(f)&locs{k}.frame<frames(f+1);
                            sum(indf)
                            fn=fieldnames(locs{k});
                            for l=1:length(fn)
                                loc2{k}.(fn{l})=locs{k}.(fn{l})(indf);
                            end
                        end
                        stats=make_statistics2(loc2,p,false);
                        for uf=1:length(usefields)
                            fh=usefields{uf};
                            if isfield(stats,fh{1})
                            stath=stats.(fh{1}).(fh{2});
                            for l=1:length(stath)
                                ps.(fh{1}).(fh{2}){l}(f)=stath(l);
                            end
                            end
                        end
                    end
                    filns=frames(1:end-1);
                    xl='frame';
            end
             %plot
            axall=obj.initaxis('all');
            delete(axall.Children);
            hold on
            if strcmp(xl,'filenumber')
               ax0=obj.initaxis('falloff');
                hold off
                plot(ax0,filns,ps.frames.falloff{1})
                plot(axall,filns,ps.frames.falloff{1}/max(ps.frames.falloff{1}))

                for k=2:length(locs)
                hold on
                plot(ax0,filns,ps.frames.falloff{k})
                plot(axall,filns,ps.frames.falloff{k}/max(ps.frames.falloff{k}))
                end
                xlabel('filenumber')
                ylabel('falloff frame')
                legend(modetxt)
            end
            
            ax1=obj.initaxis('number of localizations');
            hold off
            plot(ax1,filns,ps.photons.Nloc{1})
            plot(axall,filns,ps.photons.Nloc{1}/max(ps.photons.Nloc{1}))
            for k=2:length(locs)
            hold on
            plot(ax1,filns,ps.photons.Nloc{k})
            plot(axall,filns,ps.photons.Nloc{k}/max(ps.photons.Nloc{k}))
            end
            xlabel(xl)
            ylabel('number of localizations')
            legend(modetxt)
            
            ax2=obj.initaxis('photons (mu)');
            hold off
            plot(ax2,filns,ps.photons.mu{1})
            plot(axall,filns,ps.photons.mu{1}/max(ps.photons.mu{1}))
            for k=2:length(locs)
            hold on
            plot(ax2,filns,ps.photons.mu{k})
            plot(axall,filns,ps.photons.mu{k}/max(ps.photons.mu{k}))
            end
            xlabel(xl)
            ylabel('photons (mu)')
            legend(modetxt)

            ax2b=obj.initaxis('locprec');
            hold off
            plot(ax2b,filns,ps.locprec.max{1})
            hold on
            plot(ax2b,filns,ps.locprec.rising{1})

            plot(axall,filns,ps.locprec.rising{1}/max(ps.locprec.rising{1}))
            plot(axall,filns,ps.locprec.max{1}/max(ps.locprec.max{1}))
            for k=2:length(locs)
            hold on
            plot(ax2b,filns,ps.locprec.max{k})
            plot(ax2b,filns,ps.locprec.rising{k})
            plot(axall,filns,ps.locprec.rising{k}/max(ps.locprec.rising{k}))
            plot(axall,filns,ps.locprec.max{k}/max(ps.locprec.max{k}))                   

            end
            xlabel(xl)
            ylabel('locprec (nm)')
            
            for k=1:length(modetxt)
                m2{2*k-1}=[modetxt{k} ' max'];
                m2{2*k}=[modetxt{k} ' rising'];
            end
            legend(m2)

            ax3=obj.initaxis('on-time');
            hold off
            plot(ax3,filns,ps.lifetime.mu{1})
            plot(axall,filns,ps.lifetime.mu{1}/max(ps.lifetime.mu{1}))
            for k=2:length(locs)
            hold on
            plot(ax3,filns,ps.lifetime.mu{k})
            plot(axall,filns,ps.lifetime.mu{k}/max(ps.lifetime.mu{k}))
            end
            xlabel(xl)
            ylabel('on-time (frames) ')
            legend(modetxt)
            
             ax4=obj.initaxis('BG');
            hold off
            plot(ax4,filns,ps.background.max{1})
            plot(axall,filns,ps.background.max{1}/max(ps.background.max{1}))
            for k=2:length(locs)
            hold on
            plot(ax4,filns,ps.background.max{k})
            plot(axall,filns,ps.background.max{k}/max(ps.background.max{k}))
            end
            xlabel(xl)
            ylabel('max background')
            legend(modetxt)

            if isfield(ps,'PSFxnm')
                ax4=obj.initaxis('PSF');
                hold off
                plot(ax4,filns,ps.PSFxnm.max{1})
                for k=2:length(locs)
                hold on
                plot(ax4,filns,ps.PSFxnm.max{k})
                end
                xlabel(xl)
                ylabel('max PSFx (nm)')
                legend(modetxt)
            end
            xlabel(axall,xl)
            ylabel(axall,'all curves, normalized')
            
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end




function pard=guidef(obj)
pard.useroi.object=struct('String','use Roi','Style','checkbox','Value',1);
pard.useroi.position=[1,1];

pard.filter.object=struct('String','use layers/filters','Style','checkbox','Value',1);
pard.filter.position=[1,2];

pard.overview.object=struct('String','plot overview','Style','checkbox','Value',0);
pard.overview.position=[1,3];

pard.tphot.object=struct('String','photon range:','Style','text');
pard.tphot.position=[3,1];
pard.tphot.Width=1.5;

pard.photrange.object=struct('String','200 10000','Style','edit');
pard.photrange.position=[3,2.5];

pard.tlt.object=struct('String','lifetime range (frames):','Style','text');
pard.tlt.position=[4,1];
pard.tlt.Width=1.5;

pard.lifetimerange.object=struct('String','1 30','Style','edit');
pard.lifetimerange.position=[4,2.5];

% pard.checkphot.object=struct('String','use manual photon range','Style','checkbox','Value',0);
% pard.checkphot.position=[2,1];
% pard.checkphot.Width=2;
% 
% pard.photrange.object=struct('String','0','Style','edit');
% pard.photrange.position=[2,3];
p(1).value=1; p(1).on={}; p(1).off={'t1','framewindows'};
p(2).value=2; p(2).on=p(1).off; p(2).off={};
pard.timefield.object=struct('String',{{'files','frames'}},'Style','popupmenu','Callback',{{@obj.switchvisible,p}});
pard.timefield.position=[5,1];

pard.t1.object=struct('String','# time windows (for frame)','Style','text','Visible','off');
pard.t1.position=[6,1];
pard.t1.Width=2;
pard.framewindows.object=struct('String','10','Style','edit','Visible','off');
pard.framewindows.position=[6,3];

pard.plugininfo.name='Statistics vs frame/file';
pard.plugininfo.description=sprintf(['statsVsTime calculates localization statistics in dependence of the filenumber or of the frame'...
    'Locstatistics calculates all kind of statistics for localization data.\n'...
    'photons: N: number of localizations. <P>: mean. r: ratio between number of localizations above 2000 and between 1000 and 2000. mu: decay constant of exponential fit. \n'...
    'locprec: max, median and position of rising edge. \n'...
    'lifetime: how many frames does a fluorophore live (from grouping). mu: from exponential fit.\n'...
    'background: mean\n'...
    'either znm or PSFxnm.\n'...
    'locprecznm \n'...
    'frames: number of lcoalizations vs. frame. To see when localizations drop off.']);
pard.plugininfo.type='ProcessorPlugin';
end