classdef ExportEvaluationsTable<interfaces.DialogProcessor&interfaces.SEProcessor
    methods
        function obj=ExportEvaluationsTable(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            if nargin>0
                obj.handle=varargin{1};
            end
            obj.inputParameters={};
            updateplugin(0,0,obj)
        end
        
        function out=run(obj,p)  
          out=exportresults(obj,p);
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function out=exportresults(obj,p)
sites=obj.locData.SE.sites;
pluginname=p.pluginselection.selection;
results=fieldnames(sites(1).evaluation.(pluginname));
filenames=obj.getPar('filelist_long').String;
for k=1:length(sites)
    for l=1:length(results)
        if ~isfield(sites(k).evaluation.(pluginname),results{l})
            continue
        end
        varh=sites(k).evaluation.(pluginname).(results{l});
        goodval=ischar(varh) || numel(varh)==1;
        goodval=goodval & ~isstruct(varh);
        if goodval
            outstruc.(results{l}){k,1}=varh;
        else
            % outstruc.(results{l}){k}=[];
        end
    end
    outstruc.site{k,1}=sites(k).ID;
    outstruc.cell{k,1}=sites(k).info.cell;
    outstruc.filename{k,1}=filenames{sites(k).info.filenumber};
end

t=struct2table(outstruc);
file=obj.getPar('lastSMLFile');
outf=strrep(file,'_sml.mat',['_' pluginname '.csv']);
[file, path]=uiputfile(outf);
if file
    writetable(t,[path file])
end
out=t;
end


function updateplugin(a,b,obj)
if ~ isempty(obj.locData.SE.sites)
pluginnames=fieldnames(obj.locData.SE.sites(1).evaluation);
obj.guihandles.pluginselection.String=pluginnames;
end
end

function pard=guidef(obj)

pard.pluginselection.object=struct('String',{{"0"}},'Style','popupmenu');
pard.pluginselection.position=[1,1];
pard.pluginselection.Width=2;
pard.updateplugins.object=struct('String','update','Style','pushbutton','Callback',{{@updateplugin,obj}});
pard.updateplugins.position=[1,3];
pard.updateplugins.Width=1;

pard.plugininfo.type='ROI_Analyze';
end

