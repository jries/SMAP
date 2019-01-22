classdef locsfromSE<interfaces.DialogProcessor
    methods
        function obj=locsfromSE(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
                obj.inputParameters={'mainfile'};
        end
        
        function out=run(obj,p)
            if p.numberOfSubsets > 0
                numberOfSubsets = p.numberOfSubsets;
                histogram = cell(numberOfSubsets, 1);
            else
                numberOfSubsets = 1;
            end
            
            
            for i = 1:numberOfSubsets
                sites=obj.locData.SE.sites;
                induse=(getFieldAsVector(sites,'annotation','use'));
                listn=p.c_list;
                if ~strcmp(listn.selection,'all') %select only those with the right list entry
                    listv=getFieldAsVector(sites,'annotation',listn.selection,'value');
                    indlist=listv==p.c_listvalue;
                    induse=induse&indlist; 
                end
                
                induse=find(induse);
                if p.numberOfSites_perSub > 0
                    nIndus = length(induse);
                    indSelect = randi(nIndus, p.numberOfSites_perSub, 1);
                    induse=induse(indSelect);
                end
                
                nums=length(induse);
                locprecnm=zeros(nums,1);
                localizations=zeros(nums,1);
                psf=zeros(nums,1);
                switch p.c_groupfield.Value
                    case 1
                        gfield='Nlocs';
                    case 2
                        gfield='Nlocsg';
                    case 3
                        gfield='noblink';
                end
                ind=1;
                for k=induse
                    if ~p.c_listofcellsc|| any(sites(k).info.cell==p.c_listofcells)
                        info=sites(k).evaluation.countingStatistics;
                        locprecnm(ind)=info.locprecnm;
                        psf(ind)=info.PSF;
                        localizations(ind)=info.(gfield);  
                        ind=ind+1;
                        filenumbersave=sites(k).info.filenumber;
                    end
                end
                locprecnm(ind:end)=[];
                psf(ind:end)=[];
                localizations(ind:end)=[];
                length(psf)

                nonnan=~isnan(psf)&~isnan(locprecnm);
                indgood=psf>p.c_PSFmin&psf<p.c_PSFmax&nonnan;
                %list
             

                
                
                initaxis(p.resultstabgroup,'scatter');
                subplot(2,2,1);
                try
                dscatter(psf(nonnan),localizations(nonnan))
                hold on
                psfb=min(psf):5:max(psf);
                plot(psfb,bindata(psf(nonnan),localizations(nonnan),psfb,'median'))
                hold off
                xlabel('sigmapsf')
                ylabel('locs')


                subplot(2,2,2)
                dscatter(locprecnm(nonnan),localizations(nonnan))
                hold on
                lb=min(locprecnm):2:max(locprecnm);
                plot(lb,bindata(locprecnm(nonnan),localizations(nonnan),lb,'median'))
                hold off
                xlabel('locp')
                ylabel('locs')
                catch
                end


                initaxis(p.resultstabgroup,'histogram');
                dlocs=localizations(indgood);
                maxhist=myquantile(dlocs,0.995);
                % figure(43)
    %             dloc=round(maxhist/length(dlocs)*2);
                dloc=1;
                [y,x]=hist(dlocs,0:dloc:maxhist+dloc);
                xlim([1,100])
                plot(x,y)
                hmean=num2str(mean(dlocs(dlocs<maxhist)));
                hstd=num2str(std(dlocs(dlocs<maxhist)));
                hmedian=num2str(median(dlocs(dlocs<maxhist)));
                nlocs=num2str(sum(dlocs<maxhist));
                lastfile=obj.getPar('lastSMLFile');
                title(['mean: ' (hmean) ' ,std: ' (hstd)  ' ,median: ' (hmedian) ' ,nlocs: ' (nlocs)]);
                filen=obj.locData.SE.files(filenumbersave).name;
                clipboard('copy',[filen sprintf(['\t' lastfile '\t' hmean '\t' hstd '\t' hmedian '\t' nlocs])])
                % sum(y)
                % length(cluster)

    %             y=y/sum(y(:));
                if p.numberOfSubsets > 0
                    histogram{i}.h=y;
                    histogram{i}.c=x;
                else
                    histogram.h=y;
                    histogram.c=x;
                end
            end
            obj.setResults('counting_histogram',histogram);
            out=histogram;
    %             obj.locData.guiData.counting.histogram=histogram;
        end
        
%         function refit_callback(obj)
%         end
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
    end
end




function pard=guidef(obj)
% 
pard.text1.object=struct('String','PSFsigma (nm) min:max','Style','text');
pard.text1.position=[1,1];
pard.text1.Width=2;
pard.c_PSFmin.object=struct('String','100','Style','edit');
pard.c_PSFmin.position=[1,3];
pard.c_PSFmin.Width=0.5;
pard.c_PSFmax.object=struct('String','130','Style','edit');
pard.c_PSFmax.position=[1,3.5];
pard.c_PSFmax.Width=0.5;

pard.c_groupfield.object=struct('String',{{'ungrouped','grouped','blink remove'}},'Style','popupmenu');
pard.c_groupfield.position=[3,1];
pard.c_groupfield.Width=1.5;

pard.c_list.object=struct('String',{{'all','list1','list2','list3','list4'}},'Style','popupmenu');
pard.c_list.position=[3,3];
pard.c_list.Width=1;
pard.c_listvalue.object=struct('String','1','Style','edit');
pard.c_listvalue.position=[3,4];
pard.c_listvalue.Width=1;

pard.c_listofcellsc.object=struct('String','use only these cells:','Style','checkbox','Value',0);
pard.c_listofcellsc.position=[4,1];
pard.c_listofcellsc.Width=1.5;
pard.c_listofcells.object=struct('String','1:10','Style','edit');
pard.c_listofcells.position=[4,2.5];
pard.c_listofcells.Width=1.5;

pard.savehist.object=struct('String','Save histogram','Style','pushbutton','Callback',{{@obj.savehist_callback}});
pard.savehist.position=[5,1];


pard.t_numberOfSites_perSub.object = struct('String', '#of sites/subset', 'Style', 'text');
pard.t_numberOfSites_perSub.position = [6,1];
pard.t_numberOfSites_perSub.width = 1;

pard.numberOfSites_perSub.object = struct('String', '0', 'Style', 'edit');
pard.numberOfSites_perSub.position = [6,2];
pard.numberOfSites_perSub.width = 0.5;

pard.t_numberOfSubsets.object = struct('String', '#of subsets', 'Style', 'text');
pard.t_numberOfSubsets.position = [7,1];
pard.t_numberOfSubsets.width = 1;

pard.numberOfSubsets.object = struct('String', '0', 'Style', 'edit');
pard.numberOfSubsets.position = [7,2];
pard.numberOfSubsets.width = 0.5;

pard.plugininfo.name='histogram from SE';
pard.plugininfo.type='ProcessorPlugin';
end