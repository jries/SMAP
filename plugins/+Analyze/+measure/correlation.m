classdef correlation<interfaces.DialogProcessor
   % also look at FFT
   % empty / add button to add up FFT, AC
   % include 2D analysis? put 1D to ROI?
   % PERPL: calculate 2D g(r)? fit with lines or grid (2 distances)?
   properties
       currentcorr
       averagecorr
       currentpair
       averagepair
   end
    methods
        function obj=correlation(varargin)        
            obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'sr_layerson','linewidth_roi','znm_min','znm_max','sr_pixrec','layernames'};
            obj.showresults=true;
        end
        
        function out=run(obj,p)   
            out=[];
            roih=obj.getPar('sr_roihandle');
            if contains(p.mode.selection,"1D correlation")
                
                if ~contains(class(roih),'imline')
                    warning('ROI needs to be imline. cannot execute correlation analysis');
                else
                    axfftav=obj.initaxis('FFTav');
                    axcorrav=obj.initaxis('corr1Dav');
                    axcorr=obj.initaxis('corr1D');
                    axfft=obj.initaxis('FFT');
                    axprof=obj.initaxis('prof1D');
                    hold(axprof,'off')    
                    hold(axcorr,'off')
                    hold(axfft,'off') 
                    hold(axcorr,'off')
                    hold(axfftav,'off')
                end
            elseif contains(p.mode.selection,"2D pair correlation")
                axpcf=obj.initaxis('PCF');
                hold(axpcf,'off')     
            end
            if ~p.setbinwidth
                p.binwidth=p.sr_pixrec;
            end
            layers=find(p.sr_layerson);
            colors=lines(length(layers));
            locsl={};
            for k=1:length(layers)
                [locs,~, hroi]=obj.locData.getloc({'xnm','ynm','znm','locprecnm','locprecznm','xnmline','ynmline'},'layer',layers(k),'position','roi');
                % paircorrelationfunction(locs.xnm,locs.ynm,roihandle=obj.getPar('sr_roihandle'))
                if contains(p.mode.selection,"2D pair correlation")
                    locsl{k}=locs; %collected, correlations are calculated after the loop
                elseif contains(p.mode.selection,"1D correlation")
                    if  contains(class(roih),'imline')
                        ca(k)=correlationtools(locs,p.binwidth,periodguess=p.period,maxcorr=p.maxr);
                        ca(k).plot('profile',axis=axprof,color=colors(k,:))
                        ca(k).plot('autocorrelation',axis=axcorr,color=colors(k,:));
                        ca(k).plot('fft',axis=axfft,color=colors(k,:));
                        
                        if ~isempty(obj.averagecorr)
                            obj.averagecorr(k).plot('autocorrelationav',axis=axcorrav,color=colors(k,:),average=true);
                            obj.averagecorr(k).plot('fftav',axis=axfftav,color=colors(k,:),average=true);
                        end
                        
                    else
                        warning('line correlation needs line ROI')
                    end
                      
                end
            end

            if contains(p.mode.selection,"2D pair correlation")
                %auto-correlation of every layer and cross-correlation of every pair of layers
                layernames=obj.getPar('layernames');
                if isempty(layernames) %fallback if the layers are not named
                    layernames=arrayfun(@(x) ['layer' num2str(x)],1:max(layers),'UniformOutput',false);
                end
                haslocs=cellfun(@(x) numel(x.xnm)>=2,locsl);
                if ~any(haslocs)
                    warning('no localizations found. Draw a ROI that contains localizations.')
                    return
                elseif ~all(haslocs)
                    warning('no localizations in layer(s) %s, they are skipped.',strjoin(layernames(layers(~haslocs)),', '))
                    layers=layers(haslocs);locsl=locsl(haslocs);colors=colors(haslocs,:);
                end
                gr={};names={};cols=[];styles={};
                for k=1:length(layers)
                    for l=k:length(layers)
                        if k==l
                            [~,r,g]=paircrosscorrelationfunction(locsl{k}.xnm,locsl{k}.ynm,[],[],...
                                roihandle=roih,dr=p.binwidth,rmax=p.maxr);
                            names{end+1}=layernames{layers(k)};
                            cols(end+1,:)=colors(k,:);
                            styles{end+1}='-';
                        else
                            [g,r]=paircrosscorrelationfunction(locsl{k}.xnm,locsl{k}.ynm,...
                                locsl{l}.xnm,locsl{l}.ynm,roihandle=roih,dr=p.binwidth,rmax=p.maxr);
                            names{end+1}=[layernames{layers(k)} ' x ' layernames{layers(l)}];
                            cols(end+1,:)=(colors(k,:)+colors(l,:))/2;
                            styles{end+1}='--';
                        end
                        gr{end+1}=g;
                    end
                end
                for k=1:length(gr) %r(1)=0 contains the self-pairs of the auto-correlations
                    plot(axpcf,r(2:end),gr{k}(2:end),styles{k},'Color',cols(k,:));
                    hold(axpcf,'on')
                end
                hold(axpcf,'off')
                xlabel(axpcf,'r (nm)')
                ylabel(axpcf,'pair correlation function (norm)')
                legend(axpcf,names)
                out.r=r;out.g=gr;out.names=names;
            end

            if contains(p.mode.selection,"1D correlation")
                obj.currentcorr=ca;
            end

            if length(layers)>1 && contains(p.mode.selection,"1D correlation") %CC, only 2 layers can be active
                % obj.currentcorr=ca(1);
                ca(1).calculatecrosscorrelation(ca(2));
                ca(1).plot('crosscorrelation',axis=axcorr,color='k');
                if ~isempty(obj.averagecorr)
                    obj.averagecorr(1).plot('crosscorrelationav',axis=axcorrav,color='k',average=true);
                end
            end
        end
        function button_callback(obj,a,b)
            
            if contains(a.String,'Add')
                if isempty(obj.averagecorr)
                    obj.averagecorr=obj.currentcorr;
                end
                for k=1:length(obj.currentcorr)
                    obj.averagecorr(k).add(obj.currentcorr(k));
                end
                
            else %clear
                obj.averagecorr=[];
            end
           
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end



function pard=guidef(obj)
pard.inputParameters={'sr_layerson','linewidth_roi','sr_pixrec'};

 p(1).value=0; p(1).on={}; p(1).off={'binwidth'};
p(2).value=1; p(2).on={'binwidth'}; p(2).off={};

pard.mode.object=struct('String',{{'1D correlation','2D pair correlation'}},'Style','popupmenu');
pard.mode.position=[1,1];
pard.mode.Width=1.5;

pard.setbinwidth.object=struct('String','binwidth (nm) (empty: pixelsize):','Style','checkbox','Value',1,'Callback',{{@obj.switchvisible,p}});
pard.setbinwidth.position=[2,1];
pard.setbinwidth.Width=3;

pard.binwidth.object=struct('String','10','Style','edit');
pard.binwidth.position=[2,3];
pard.binwidth.Width=0.5;
pard.binwidth.TooltipString='Binwidth for profiles. If not checked, use pixel size of reconstruction';
pard.setbinwidth.TooltipString=pard.binwidth.TooltipString;


pard.periodt.object=struct('String','approximate period (nm)','Style','text');
pard.periodt.position=[3,1];
pard.periodt.Width=1.5;
pard.period.object=struct('String','200','Style','edit');
pard.period.position=[3,2.5];
pard.period.Width=0.5;

pard.maxrt.object=struct('String','max correlation (nm)','Style','text');
pard.maxrt.position=[3,3];
pard.maxrt.Width=1.5;
pard.maxr.object=struct('String','1000','Style','edit');
pard.maxr.position=[3,4.5];
pard.maxr.Width=0.5;

pard.add.object=struct('String','Add','Style','pushbutton','Callback',@obj.button_callback);
pard.add.position=[4,1];
pard.add.Width=0.5;

pard.clear.object=struct('String','Clear','Style','pushbutton','Callback',@obj.button_callback);
pard.clear.position=[4,4];
pard.clear.Width=0.5;

% pard.text2.object=struct('String','fitmodel:','Style','text');
% pard.text2.position=[2,1];
% 
% pard.fitmodel.object=struct('String','Gauss|Flat|Disk|Ring|Distance','Style','popupmenu');
% pard.fitmodel.position=[2,2];
% 
% pard.restrictsigma.object=struct('String','sigma=<locp>','Style','checkbox');
% pard.restrictsigma.position=[2,3];
% 
% 
% p(1).value=0; p(1).on={}; p(1).off={'linelength'};
% p(2).value=1; p(2).on={'linelength'}; p(2).off={};
% pard.linelengthcheck.object=struct('String','set length (nm)','Style','checkbox','Callback',{{@obj.switchvisible,p}});
% pard.linelengthcheck.position=[3,1];
% 
% pard.linelength.object=struct('String','250','Style','edit');
% pard.linelength.position=[3,2];
% pard.linelength.TooltipString='This overrides the length of the ROI and uses a well-defined ROI. Useful for direct comparison.';
% pard.linelengthcheck.TooltipString=pard.linelength.TooltipString;
pard.plugininfo.name='Correlation';
pard.plugininfo.description=sprintf('correlation analysis');
pard.plugininfo.type='ProcessorPlugin';
end