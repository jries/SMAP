classdef BatchAnalysis<interfaces.DialogProcessor
% Batch analysis of many files with defined plugins. Create a Batchanlaysis
% Tab in the Anlaysis tab and add all plugins for batch analysis there.';

    methods
        function obj=BatchAnalysis(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.showresults=false;
        end
        
        function out=run(obj,p)         
            out=[];
            gsmap=obj.getPar('mainGui');
            allpars=gsmap.getGuiParameters(true);
            
            switch p.loadwhat.Value
                case 1 %txt file with file list
                    infile=p.filelistfile;
                    fID=fopen(infile);
                    fh=fgetl(fID);
                    k=1;
                    while fh>0
                        filelist{k}=fh;
                        k=k+1;
                         fh=fgetl(fID);
                    end
                    fclose(fID);
                case 2 %directory
                    files = dir([p.filelistfile filesep '**' filesep '*_sml.mat']);
                    for k=1:length(files)
                        filelist{k}=[files(k).folder filesep files(k).name];
                    end
            end
            gabatch=gsmap.children.guiAnalysis.children.Batchanalysis;
            aplugins=fieldnames(gabatch.children);
            gFile=gsmap.children.guiFile;
            gFormat=gsmap.children.guiRender.children.guiFormat;
            
            [outp,outf]=fileparts(p.outdir);
            results=[];
            
            
            for a=1:length(aplugins)
                mkdir([outp filesep outf filesep aplugins{a}]);
            end
            for f=1:length(filelist)
                [path,file,ext]=fileparts(filelist{f});
                gFile.loadbutton_callback( 0,0,0,[path filesep],[file ext]);
                gsmap.setGuiParameters(allpars,true);
                gFormat.resetview;
                %re-evaluate if needed
                if p.evalutesites
                    gsmap.children.guiSites.children.eval.redrawall;
                end
                %here set FoV, or reset. Maybe set ROI
                for a=1:length(aplugins)
                   ahere=gabatch.children.(aplugins{a});
%                    try
                       re=ahere.processgo;
                       re
                       results.(aplugins{a}){f}=re;
                       outfig=ahere.resultstabgroup.Parent;
                       savefig(outfig,[outp filesep outf filesep aplugins{a} filesep file '.fig']);
%                    catch err
%                        disp('output could not be saved')
%                        err
%                    end
                end
            end
            save(p.outdir,'results','filelist')

% - Text file with all SML files
% if this is directory (maybe with switch): recursive search
% - Make new tab (batch) in analyse, add all analysis plugins, configure
% - Load first file, filter, set everything
% - Evaluation plugins: configure.
% 
% SMAP:
% - Read all GUI parameters
% - Load file without parameters
% - Restore all GUI parameters
% - If needed: evaluate all
% - Loop over all analysis plugins and execute
%     - Make sure to capture outputs, also save output figures if show is active
%     - Also capture output numbers, write in csv file

        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function load_filelist(a,b,obj)
loadwhat=obj.getSingleGuiParameter('loadwhat');
f=obj.getSingleGuiParameter('filelistfile');

switch loadwhat.Value
    case 1
        if isempty(f)
            f=([fileparts('lastSMLFile'),'*.txt']);
        end
        if ~contains(f,'.txt')
            f=[f filesep '*.txt'];
        end
        [f,p]=uigetfile(f);
    case 2
        if isempty(f)
            f=([fileparts('lastSMLFile')]);
        end
        f=uigetdir(f);
        p='';
end
if f
obj.setGuiParameters(struct('filelistfile',[p f]));
end
end

function set_outdir(a,b,obj)
f=obj.getSingleGuiParameter('outdir');
if isempty(f)
    flist=obj.getSingleGuiParameter('filelistfile');
    if contains(flist,'.txt')
        f=strrep(flist,'.txt','_results.mat');
    else
        f=[flist '_results.mat'];
    end
end
[f,p]=uiputfile(f);
if f
obj.setGuiParameters(struct('outdir',[p f]));
end
end

function pard=guidef(obj)

pard.text.object=struct('Style','text','String','Create a Batchanlaysis tab in the Analysis tab and add all plugins that you want to evaluate.');
pard.text.position=[1,1];
pard.text.Width=4;

pard.loadwhat.object=struct('Style','popupmenu','String',{{'file list', 'directory'}});
pard.loadwhat.position=[2,1];
pard.loadwhat.Width=1;

pard.filelistfile.object=struct('Style','edit','String','');
pard.filelistfile.position=[2,2];
pard.filelistfile.Width=2.5;

pard.loadfilelist.object=struct('Style','pushbutton','String','load','Callback',{{@load_filelist,obj}});
pard.loadfilelist.position=[2,4.5];
pard.loadfilelist.Width=0.5;



pard.outdirt.object=struct('Style','text','String','results directory');
pard.outdirt.position=[3,1];
pard.outdirt.Width=1;

pard.outdir.object=struct('Style','edit','String','');
pard.outdir.position=[3,2];
pard.outdir.Width=2.5;

pard.setoutfile.object=struct('Style','pushbutton','String','set','Callback',{{@set_outdir,obj}});
pard.setoutfile.position=[3,4.5];
pard.setoutfile.Width=0.5;

pard.evalutesites.object=struct('Style','checkbox','String','evaluate ROIs');
pard.evalutesites.position=[4,1];
pard.evalutesites.Width=1;

pard.plugininfo.type='ProcessorPlugin';

pard.plugininfo.description='Batch analysis of many files with defined plugins. Create a Batchanlaysis Tab in the Anlaysis tab and add all plugins for batch analysis there.';
end