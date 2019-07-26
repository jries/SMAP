classdef GuiFilterTable< interfaces.LayerInterface
%     Filter on any attribute of the localizations
    properties
        filter
        excludefields;
    end
    methods
        function obj=GuiFilterTable(varargin)
            obj@interfaces.LayerInterface(varargin{:});
            obj.outputParameters={'filter'};
            obj.propertiesToSave={'filter'};
        end
        function makeGui(obj)
           obj.excludefields={'peakfindxnm','peakfindynm','groupindex','photerr','inframes'};
           columnname = {'Name','min set','min','mean','max','max set','filt','inv'};
           columnformat = {'char','short g','short g','short g','short g','short g','logical','logical'};
           columnedit=[false, true, false, false, false, true, true, true];
           
           s{1,8}=false;
           h.table = uitable('Parent',obj.handle,'Data', s,... 
            'ColumnName', columnname,'Units','normalized',...
            'ColumnFormat', columnformat,...
            'ColumnEditable', columnedit,'FontSize',12,...
            'RowName',[],'CellSelectionCallback',{@cellSelect_callback,obj},'CellEditCallback',{@cellEdit_callback,obj});
           h.table.Position(3) = .95;
           h.table.Position(4) = .95;
           h.table.Position(2) = 0.025;
           h.table.Units='pixel';
           ppix=h.table.Position;
           h.table.Units='normalized';
           w=ppix(3)*.95;
           widths=[2 1 1 1 1 1 0.5 0.5];
           widths=widths/sum(widths);
           h.table.ColumnWidth=num2cell(w*widths);
           h.handle=obj.handle;
           obj.guihandles=h;
           obj.initGui;
           callobj=obj;
           obj.addSynchronization('locFields',[],[],{@callobj.updateGui})          
           inittable={{'xnm',[-inf inf],true},{'ynm',[-inf inf],true},...
               {'filenumber',[-inf inf],true},{'locprecnm',[0 31],true}};
           for k=1:length(inittable)
                obj.filter.(inittable{k}{1}).minmax=inittable{k}{2};
                obj.filter.(inittable{k}{1}).auto=inittable{k}{3};
           end        
           obj.addSynchronization([obj.layerprefix 'selectedField'],[],[],{@callobj.selectedField_callback})           
           if ~isempty(obj.locData.loc)
            updateGui(obj,0,0)
           end
        end
        
        function setGuiParameters(obj,p,setchildren)
            setGuiParameters@interfaces.GuiModuleInterface(obj,p,setchildren);  
            obj.guihandles.table.Data=p.table.Data;
            obj.setPar('filtertable',p.table.Data,'layer',obj.layer);
            refilter(obj);
        end

        function selectedField_callback(obj)
            
            sfield=obj.getPar('selectedField','layer',obj.layer);
            if isempty(sfield)
                return
            end
            field=sfield{1};

            fmin=sfield{2};
            fmax=sfield{3};
            obj.filter.(field).minmax=[fmin fmax];
            if length(sfield)<4||isempty(sfield{4})
                fauto=[];
            else
                fauto=sfield{4}(1);
                obj.filter.(field).auto=fauto;
            end
            s=obj.getPar('filtertable','layer',obj.layer);
            sold=s;            

            if ~isempty(s)
            indf=find(strcmp(s(:,1),field));
            else
                indf=[];
            end
            
            if ~isempty(indf)
                if length(sfield)>3&&length(sfield{4})>1
                    s{indf,8}=logical(sfield{4}(2));
                    obj.filter.(field).invert=s{indf,8};
                else
                    
                end
            
                s{indf,2}=fmin;
                s{indf,6}=fmax;
                
                
                if ~isempty(fauto)
                    s{indf,7}=logical(fauto);
                    if ~fauto
                        obj.locData.removeFilter(field,obj.layer);
                    end
                end

                obj.guihandles.table.Data=s;
                obj.setPar('filtertable',s,'layer',obj.layer);
                if  s{indf,2}~=sold{indf,2} || s{indf,6}~=sold{indf,6}||s{indf,7}~=sold{indf,7}||s{indf,8}~=sold{indf,8}
                    refilter(obj,field)
                end
            end          
        end
        function updateGui(obj,object,event)
           obj.status('filter locs')
           locs=obj.locData.loc;
           if ~isempty(locs)
               fnall=fieldnames(locs);
               for k=1:length(fnall)
                   if ~isnumeric(locs.(fnall{k})) && ~islogical(locs.(fnall{k})) %filtering only possible on numeric vectors
                       obj.excludefields=union(obj.excludefields,fnall{k});
                   end
               end
               fn=setdiff(fnall,obj.excludefields);
               for k=1:length(fn)
                   if ~isfield(obj.filter,fn{k})||~isfield(obj.filter.(fn{k}),'auto')||~isfield(obj.filter.(fn{k}),'minmax')||~isfield(obj.filter.(fn{k}),'invert')||...
                           isempty(obj.filter.(fn{k}).auto)||isempty(obj.filter.(fn{k}).minmax)
%                    if isfield(obj.filter,fn{k})&&~isempty(obj.filter.(fn{k}).minmax)&&~isempty(obj.filter.(fn{k}).auto)
%                    else
                       obj.filter.(fn{k}).minmax=[-inf inf];
                       obj.filter.(fn{k}).auto=false;
                       obj.filter.(fn{k}).invert=false;
                   end
               end
               sss=obj.guihandles.table.Data;
               noupdate=size(sss,1)==length(fn) && all(strcmp(sss(:,1),fn));
               
               for k=1:length(fn)
                   val=obj.locData.loc.(fn{k});
                   s{k,1}=fn{k};
                   s{k,2}=obj.filter.(fn{k}).minmax(1);
                   s{k,6}=obj.filter.(fn{k}).minmax(2);
                   if ~noupdate || length(val)<2e6%this can take time, only update if fields changed.
                       s{k,3}=min(val);
                       s{k,4}=mean(val);
                       s{k,5}=max(val);
                   end
                   s{k,7}=logical(obj.filter.(fn{k}).auto);
                   s{k,8}=obj.filter.(fn{k}).invert;
               end
               obj.guihandles.table.Data=s;
                obj.setPar('filtertable',s,'layer',obj.layer);
                %apply filter to all fields not auto   
                obj.locData.layer(obj.layer).filter=[];
                obj.locData.layer(obj.layer).groupfilter=[];
                
%                 allfields={'PSFxnm','znm','locprecznm','locprecnm','frame'};
%                 for k=1:length(fn)
%                     if s{k,7}&&   any(strcmp(allfields,s{k,1}))%filter
%                         sf={s{k,1},s{k,2},s{k,6},s{k,7},true};
%                         obj.setPar('selectedField',sf,'layer',obj.layer)
%                     end
%                 end


%                for k=1:length(fn)
%                    if s{k,7}
%                        obj.locData.filter(fn{k},obj.layer);
%                    end
%                end 
                obj.locData.filter([],obj.layer);
           end
            obj.status('filter locs done')
        end
    end
end

% function ind=findrow(cell,pattern)
% ind=[];
% for k=1:length(cell)
%     if strcmp(cell{k},pattern)
%         ind=k;
%     end
% end
% end

function cellSelect_callback(object,data,obj)
if ~isempty(data.Indices)&&data.Indices(2)==1
    cellEdit_callback(object,data,obj)
end
end

function cellEdit_callback(object,data,obj)
    ind=data.Indices;
    s=obj.guihandles.table.Data;
    if ~isempty(ind)
        fn=s{data.Indices(1),1};
        sf={fn,s{data.Indices(1),2},s{data.Indices(1),6},[s{data.Indices(1),7},s{data.Indices(1),8}],true,true};
        obj.setPar('selectedField',sf,'layer',obj.layer)
        obj.selectedField_callback
    end
end

function refilter(obj,fn)
if nargin>1
    obj.locData.filter(fn,obj.layer)
else
    obj.locData.filter([],obj.layer)
end
end