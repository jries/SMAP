classdef PlotRingCME<interfaces.DialogProcessor&interfaces.SEProcessor
    properties

        results
        resultsfigure
    end
    methods
        function obj=PlotRingCME(varargin)        
            obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={};
        end
        
        function out=run(obj,p)  
            totalresults=obj.results(1);
            for k=2:length(obj.results)
                totalresults=addresults(totalresults,obj.results(k));
            end
            obj.resultsfigure=cmeRingplotresults(p,totalresults);
            out=totalresults;
%             out=obj.results;
%           analyzeRingsCME(obj.SE,p)
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function load_callback(obj,a,b,mode)
            switch mode
                case {'add'}
                    [f,path]=uigetfile('*.mat','results file');
                    if path
                        ls=load([path f]);
                        if isfield(ls,'results')
%                             disp('adding not implemented')
                            if isempty(obj.results)
                                obj.results=ls.results;
                            else
                                obj.results(end+1)=obj.results;
                            end
                            obj.results(end).filename=[path f];
                        end
                    end
                
                case 'remove'
                    p=obj.getAllParameters;
                    v=p.resultslist.Value;
                    obj.results(v)=[];
                    

            end
            setresultslist(obj);
        end
         function save_callback(obj,a,b)
            p=obj.getAllParameters;
            selection=p.savemenu.selection;
%             [~,~,ext]=fileparts(selection);
            [file,path]=uiputfile(selection);
            if path
                switch selection
                    case 'results.mat'
                        results=obj.results;
                        save([path file],'results');
                    case 'sites.mat'
                        sites=obj.sites;
                        save([path file],'sites');
                    case 'results.pdf'
                        export_fig([path file],'-pdf','-nocrop',obj.resultsfigure)
%                         saveas(obj.resultsfigure,[path file],'epsc')

                    case 'average.tif'
                        image=obj.results.sumimage/obj.results.numsites;
                        saveastiff(uint16(image/max(image(:))*2^16),[path file])
                    case 'rad_av.mat'
                        results.sumimage=obj.results.sumimage;
                        results.sumrdensity=obj.results.sumrdensity;
                        results.numberofsites=obj.results.numsites;
                        save([path file],'results');
                end
            end
         end
    end
end




function results=addresults(results,r2)
lc=length(r2.N);
    results.N(end+1:end+lc)=r2.N;
    results.rc(end+1:end+lc)=r2.rc;
    results.dr(end+1:end+lc)=r2.dr;
    results.ro(end+1:end+lc)=r2.ro;
    results.sigma(end+1:end+lc)=r2.sigma;
    maxfn=max(results.filenumber);
    results.filenumber(end+1:end+lc)=r2.filenumber+maxfn;
    results.Nnormmean(end+1:end+lc)=r2.Nnormmean;
    results.Nnormmedian(end+1:end+lc)=r2.Nnormmedian;
    
    results.rdensity(end+1:end+lc)=r2.rdensity;
    results.images(end+1:end+lc)=r2.images;
    results.ac(end+1:end+lc)=r2.ac;
    results.rdensityn(end+1:end+lc)=r2.rdensityn;
    results.acthetan(end+1:end+lc)=r2.acthetan;
    results.sitenames(end+1:end+lc)=r2.sitenames;
    
    lf=length(r2.filenumberrange);
    results.filenumberrange(end+1:end+lf)=r2.filenumberrange+maxfn;
    results.Nfilemean(end+1:end+lf)=r2.Nfilemean;
    results.Nfilemedian(end+1:end+lf)=r2.Nfilemedian;
    results.sumimage=results.sumimage+r2.sumimage;
    maxind=max(length(results.sumrdensity),length(r2.sumrdensity));
    results.sumrdensity(maxind)=0;r2.sumrdensity(maxind)=0;
    results.sumrdensity=results.sumrdensity+r2.sumrdensity;
    results.numsites=results.numsites+lc;

end
function setresultslist(obj)
filelist={obj.results(:).filename};
obj.guihandles.resultslist.String=filelist;
obj.guihandles.resultslist.Value=min(max(obj.guihandles.resultslist.Value,length(filelist)),1);
end
function pard=guidef(obj)
pard.resultslist.object=struct('String','empty','Style','listbox');
pard.resultslist.position=[5,1];
pard.resultslist.Height=5;
pard.resultslist.Width=4;

pard.add.object=struct('String','Add','Style','pushbutton','Callback',{{@obj.load_callback,'add'}});
pard.add.position=[6,1];

pard.remove.object=struct('String','Remove','Style','pushbutton','Callback',{{@obj.load_callback,'remove'}});
pard.remove.position=[6,2];

pard.savebutton.object=struct('String','Save','Style','pushbutton','Callback',{{@obj.save_callback}});
pard.savebutton.position=[6,4];

pard.savemenu.object=struct('String',{{'results.pdf','average.tif','results.mat','sites.mat','rad_av.mat'}},'Style','popupmenu');
pard.savemenu.position=[7,4];

pard.t1.object=struct('String','max dr/ro','Style','text');
pard.t1.position=[7,1];
pard.maxdrro.object=struct('String','1.5','Style','edit');
pard.maxdrro.position=[7,2];
pard.plugininfo.type='ROI_Analyze';
    
end