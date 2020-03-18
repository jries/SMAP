classdef export_coordinates<interfaces.DialogProcessor
%     exports localization data in .txt, .csv, .xls format. The user can
%     choose which localization attributes to export and if to export all
%     localizations or filtered localizations.
    properties
        exportfields={'xnm';'ynm';'znm';'frame';'locprecnm';'phot';'bg';'layer'};
    end
    methods
        
        function obj=export_coordinates(varargin)  
                obj@interfaces.DialogProcessor(varargin{:}) ;
                obj.inputParameters={'filelist_long','mainfile','mainGui','numberOfLayers','sr_layerson'};
                
        end
        
        function out=save(obj,p,defaultfilename)
            if nargin < 3
                defaultfilename=false;
            end
            
            obj.status('export localizations')
            par=p;
            
            switch par.savevisible.selection
                case 'all'
                    if par.grouped
                        gr='grouped';
                    else
                        gr='ungrouped';
                    end
                    locs=obj.locData.getloc(obj.exportfields,'grouping',gr,'position','all');
                    taball=struct2table(locs);
                otherwise
                    taball=[];
                    for k=1:p.numberOfLayers
                        if p.sr_layerson(k)
                            locs=obj.locData.getloc(obj.exportfields,'layer',k,'position','roi');
                            
                            locs.layer=k*ones(size(locs.(obj.exportfields{1})));
                            ltab=struct2table(locs);
                            taball=vertcat(taball,ltab);
                        end
                    end
            end
            
            
%             fn=p.filelist_long.selection;
            fn=obj.getPar('lastSMLFile');
            [path,file,ext]=fileparts(fn);
            of=[path filesep file  '.' par.format.selection];
            
            if ~defaultfilename
                [f,path]=uiputfile(of);
                if ~f
                    out=0;
                    return
                end
                of=[path f];
            end
            out=of;
            writetable(taball,of);
            obj.status('save done')
        end
        function initGui(obj)
            obj.guihandles.outputfields.String=sprintf('%s, ',obj.exportfields{:});
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function out=run(obj,p)
            out=obj.save(p,true);
        end
    end
end

function selectfields(a,b,obj)
p=obj.getGuiParameters;

fieldsh=fieldnames(obj.locData.loc);
if contains(p.savevisible.selection,'visible')
    fieldsh{end+1}='layer';
    obj.exportfields{end+1}='layer';
end

ind=contains(obj.exportfields,fieldsh);
fields=vertcat(obj.exportfields(ind),setdiff(fieldsh,obj.exportfields));

position=1:length(fields);
[~,ind]=intersect(fields,obj.exportfields);
exportthis=false(length(fields),1);
exportthis(ind)=true;

ds.export=exportthis;
ds.fields=fields;
ds.position=position';

data=table2cell(struct2table(ds));


h=figure('Name','Select Fields','ToolBar','none','MenuBar','none');
h.Position(3)=h.Position(3)*.6;
ht=uitable(h);
ht.Units='normalized';
ht.Position=[0 0.1 1 .9];
ht.Data=data;
ht.ColumnEditable=[true false true];
ht.ColumnName={'save','field','position'};
hb=uicontrol('Position',[5 5 50 30],'String','OK','Callback',@closef);

    function closef(a,b)
        data=ht.Data;
        close(h);
        expf=data([data{:,1}],2);
        pos=data([data{:,1}],3);
        [~,inds]=sort([pos{:}]);
        
        obj.exportfields=expf(inds);
%         expf(inds)
        obj.guihandles.outputfields.String=sprintf('%s, ',obj.exportfields{:});
    end

%layer: more than one layer: do it layer-wise. Export field: layer
end



function pard=guidef(obj)
% pard.plugininfo={'csv saver'};
pard.savevisible.object=struct('Style','popupmenu','Visible','on','String',{{'all','visible ROI'}},'Value',1);
pard.savevisible.position=[1,1];
pard.savevisible.Width=1;

pard.grouped.object=struct('Style','checkbox','Visible','on','String','grouped (if all)','Value',0);
pard.grouped.position=[1,2];
pard.grouped.Width=2;

pard.format.object=struct('Style','popupmenu','Visible','on','String',{{'csv','txt','dat','xls'}});
pard.format.position=[1,4];
pard.format.Width=1;

pard.outputfields.object=struct('Style','text','Visible','on','String','xnm, ynm, frame');
pard.outputfields.position=[2,1];
pard.outputfields.Width=3;

pard.selectoutputfields.object=struct('Style','pushbutton','Visible','on','String','select','Callback',{{@selectfields,obj}});
pard.selectoutputfields.position=[2,4];
pard.selectoutputfields.Width=1;

% pard.defaultfilename.object=struct('Style','checkbox','Visible','on','String','default filename','Value',0);
% pard.defaultfilename.position=[3,1];
% pard.defaultfilename.Width=2;

pard.plugininfo.type='SaverPlugin';
pard.plugininfo.description='exports localization data in .txt, .csv, .xls format. The user can choose which localization attributes to export and if to export all localizations or filtered localizations.';
end