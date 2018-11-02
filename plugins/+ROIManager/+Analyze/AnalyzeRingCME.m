classdef AnalyzeRingCME<interfaces.DialogProcessor&interfaces.SEProcessor
    properties
        sites
        results
        resultsfigure
    end
    methods
        function obj=AnalyzeRingCME(varargin)        
            obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={};
        end
        
        function out=run(obj,p)  
            if p.restrictcheck
                if length(p.usesites)==1
                    range=1:min(p.usesites,length(obj.SE.sites));
                else
                    range=p.usesites;
                end
            else 
                range=1:length(obj.SE.sites);
            end
            obj.sites=obj.SE.sites(range);
                
            if isempty(obj.sites)
                disp('no sites loaded')
                return
            end
            obj.results=cmeresults(obj.sites);   
            obj.resultsfigure=cmeRingplotresults(p,obj.results);
            out=obj.results;

        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end

        function save_callback(obj,a,b)
            p=obj.getAllParameters;
            selection=p.savemenu.selection;
            pf=selection;
            if strcmp(selection,'all (not sites)')
                pf='all.mat';
            end
%             [~,~,ext]=fileparts(selection);
            [file,path]=uiputfile(pf);
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
                        if size(image,3)==3
                            options.color=true;
                        else
                            options.color=false;
                        end
                        saveastiff(uint16(image/max(image(:))*2^16),[path filen '.tif'],options)
                        saveastiff(uint16(image/max(image(:))*2^16),[path file])
                    case 'rad_av.mat'
                        results.sumimage=obj.results.sumimage;
                        results.sumrdensity1=obj.results.sumrdensity1;
                        results.sumrdensity2=obj.results.sumrdensity2;
                        results.numberofsites=obj.results.numsites;
                        save([path file],'results');
                    case 'all (not sites)'
                        [~,filen]=fileparts(file);
                        obj.saveall(path,filen);
%                         results=obj.results;
%                         save([path filen '_results.mat'],'results');
%                         results=[];
%                         results.sumimage=obj.results.sumimage;
%                         results.sumrdensity1=obj.results.sumrdensity1;
%                         results.sumrdensity2=obj.results.sumrdensity2;
%                         results.numberofsites=obj.results.numsites;
%                         save([path filen '_rad_av.mat'],'results');
%                         
%                         image=obj.results.sumimage/obj.results.numsites;
%                         
%                         if size(image,3)==3
%                             options.color=true;
%                         else
%                             options.color=false;
%                         end
%                         saveastiff(uint16(image/max(image(:))*2^16),[path filen '.tif'],options)
%                         export_fig([path filen '.pdf'],'-pdf','-nocrop',obj.resultsfigure)

                end
            end
        end
        
        function saveall(obj,path,filen)
            results=obj.results;
            save([path filen '_results.mat'],'results');
            results=[];
            results.sumimage=obj.results.sumimage;
            results.sumrdensity1=obj.results.sumrdensity1;
            results.sumrdensity2=obj.results.sumrdensity2;
            results.numberofsites=obj.results.numsites;
            save([path filen '_rad_av.mat'],'results');

            image=obj.results.sumimage/obj.results.numsites;

            if size(image,3)==3
                options.color=true;
            else
                options.color=false;
            end
            saveastiff(uint16(image/max(image(:))*2^16),[path filen '.tif'],options)
            export_fig([path filen '.pdf'],'-pdf','-nocrop',obj.resultsfigure)
            
        end
    end
end

function results=cmeresults(sites)
if isfield(sites(1).evaluation.CME2DRing.circfit,'Ncirc1')
    name.Ncirc='Ncirc1';
    name.rc='r1';
    name.dr='dr1';
    name.ro='r1';
    name.sigma='sigma1';
else
    name.Ncirc='Ncirc';
    name.rc='r2D';
    name.dr='dr';
    name.ro='r2D';
    name.sigma='sigma';
end

circfitfields={'evaluation','CME2DRing','circfit'};
imfitfields={'evaluation','CME2DRing','imfit'};
results.N=getFieldAsVector(sites,circfitfields{:},name.Ncirc);
results.rc=getFieldAsVector(sites,circfitfields{:},name.rc);

results.images=getFieldAsVector(sites,imfitfields{:},'image');
results.dr=getFieldAsVector(sites,imfitfields{:},name.dr);
results.ro=getFieldAsVector(sites,imfitfields{:},name.ro);

results.sigma=getFieldAsVector(sites,imfitfields{:},name.sigma);


results.ac=getFieldAsVector(sites,imfitfields{:},'profiles1','thetaAC');
if isfield(sites(1).evaluation.CME2DRing.circfit,'profiles2')
    results.rdensity1=getFieldAsVector(sites,imfitfields{:},'profiles1','rdensity');
    results.rdensity2=getFieldAsVector(sites,imfitfields{:},'profiles2','rdensity');
    rall=cell2mat(results.rdensity1');
    stdr=std(rall,1);
    results.rdensitystd1=stdr;
    rall=cell2mat(results.rdensity2');
    stdr=std(rall,1);
    results.rdensitystd2=stdr;
else
    results.rdensity1=getFieldAsVector(sites,imfitfields{:},'profiles1','rdensity');
    rall=cell2mat(results.rdensity1');
    stdr=std(rall,1);
    results.rdensitystd1=stdr;
end
results.rdensityn=getFieldAsVector(sites,imfitfields{:},'profiles1','rn');
results.acthetan=getFieldAsVector(sites,imfitfields{:},'profiles1','thetan');

results.filenumber=getFieldAsVector(sites,'info','filenumber');
results.sitenames=getFieldAsVector(sites,'name');
    
    
filenumberrange=1:max(results.filenumber);
results.filenumberrange=filenumberrange;
Nfilemean=0;
Nfilemedian=0;
% Nnorm=N;
N=results.N;
Nmed=median(N);
Nmean=mean(N);
for k=filenumberrange
    ink=results.filenumber==k;
    Nfilemean(k)=mean(N(ink));
    Nfilemedian(k)=median(N(ink));
    Nnormmean(ink)=N(ink)/Nfilemean(k)*Nmean;
    Nnormmedian(ink)=N(ink)/Nfilemedian(k)*Nmed;  
end

results.Nfilemean=Nfilemean;
results.Nfilemedian=Nfilemedian;
results.Nnormmean=Nnormmean;
results.Nnormmedian=Nnormmedian;

results.numsites=length(sites);
sumim=zeros(size(results.images{1}));
sumrdensity1=zeros(size(results.rdensity1{1}));

if ~isfield (results, 'rdensity2')
    sumrdensity2=sumrdensity1;
else
    sumrdensity2=zeros(size(results.rdensity2{1}));
end

for k=1:length(sites)
%     size(results.images{k})
    sumim=results.images{k}+sumim;
    sumrdensity1=results.rdensity1{k}+sumrdensity1;
    if ~isfield (results, 'rdensity2')
        sumrdensity2=sumrdensity1;
    else
        sumrdensity2=results.rdensity2{k}+sumrdensity2;
    end
end
results.sumimage=sumim;
results.sumrdensity1=sumrdensity1;
results.sumrdensity2=sumrdensity2;
results.sumrdensityn=results.rdensityn{1};


% add intensityTiff eval

if isfield(sites(1).evaluation,'intensityTiff')
    results.intensityTiff=getFieldAsVector(sites,'evaluation','intensityTiff');
end

end


function pard=guidef(obj)


pard.savebutton.object=struct('String','Save','Style','pushbutton','Callback',{{@obj.save_callback}});
pard.savebutton.position=[1,1];

pard.savemenu.object=struct('String',{{'all (not sites)','results.pdf','average.tif','results.mat','sites.mat','rad_av.mat'}},'Style','popupmenu');
pard.savemenu.position=[2,1];

pard.t1.object=struct('String','max dr/ro','Style','text');
pard.t1.position=[1,3];
pard.maxdrro.object=struct('String','1.5','Style','edit');
pard.maxdrro.position=[1,4];

pard.restrictcheck.object=struct('String','only sites','Style','checkbox','Value',0);
pard.restrictcheck.position=[2,3];
pard.usesites.object=struct('String','150','Style','edit');
pard.usesites.position=[2,4];
pard.usesites.TooltipString='first N sites, or LIst of sites or s1:s2 notation';
pard.plugininfo.type='ROI_Analyze';
end