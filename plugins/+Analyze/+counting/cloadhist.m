classdef cloadhist<interfaces.DialogProcessor
    methods
        function obj=cloadhist(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
                obj.inputParameters={'mainfile'};
        end
        
        function out=run(obj,p)
%             hax=initaxis(p.resultstabgroup,'histograms');
            hthere=obj.getResults('counting_histogram');
            f=figure(90);
            clf(f)
            hax=gca;
            hax.FontSize=16;
            hax.XLabel.String='Number of localizations';
            hax.YLabel.String='Frequency';
            hax.Title.String='Averaged histogram of localizations';
            plot(hthere.c,hthere.h/sum(hthere.h),'g','Parent',hax);
            hold on
            
            of=p.mainfile;
            if isempty(of)
                off='*_hist.mat';
            else
                [path,f,ext]=fileparts(of);
                off=[path filesep f '_hist.mat'];
            end
            [file,path]=uigetfile(off,'MultiSelect','on');
            if path
                if ~iscell(file)
                    file={file};
                end
                for k=1:length(file)
                    l=load([path file{k}]);
                    hload=l.histogram;
                     plot(hload.c,hload.h/sum(hload.h),'Parent',hax);
                    if length(hload.c)<length(hthere.c)
                        hthere.h(1:length(hload.c))= hthere.h(1:length(hload.c))+hload.h;
                    else
                        hloadh=hload.h(1:length(hthere.c));
                        hload.h(1:length(hthere.c))=hloadh(:)+hthere.h(:);
                        hthere=hload;
                    end
                end
            end
            plot(hthere.c,hthere.h/sum(hthere.h),'k','Parent',hax,'LineWidth',2);
            hold off
            obj.setResults('counting_histogram',hthere);
            out=hthere;

        end
        
        function refit_callback(obj)
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function savehist_callback(obj,a,b)
            p=obj.getAllParameters;
            histogram=obj.getResults('counting_histogram');
            of=p.mainfile;
            if isempty(of)
                off='*_hist.mat';
            else
                [path,f,ext]=fileparts(of);
                off=[path filesep f '_hist.mat'];
            end
            [f,path]=uiputfile(off);
            if f
                save([path,f],'histogram')
            end
            
        end
        function empty_callback(obj,a,b)
            h.c=[];h.h=[];
            obj.setResults('counting_histogram',h);
        end
    end
end




function pard=guidef(obj)

pard.emptyhist.object=struct('String','empty histogram','Style','pushbutton','Callback',{{@obj.empty_callback}});
pard.emptyhist.position=[2,1];

pard.savehist.object=struct('String','Save histogram','Style','pushbutton','Callback',{{@obj.savehist_callback}});
pard.savehist.position=[5,1];


pard.plugininfo.name='load brightness histograms';
pard.plugininfo.type='ProcessorPlugin';
end