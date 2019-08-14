classdef GroupAny<interfaces.DialogProcessor
%     Combines localizations which have the same value in a selected field
    methods
        function obj=GroupAny(varargin)      
            obj@interfaces.DialogProcessor(varargin{:}) ;  
        end
        
        function out=run(obj,p)
            out=[];
            obj.setPar('undoModule','GroupAny');
            notify(obj.P,'backup4undo');
            gfield=p.locfield.selection;
            gval=obj.locData.loc.(gfield);
            sortmatrix=horzcat(single(obj.locData.loc.filenumber),obj.locData.loc.channel);
            if isfield(obj.locData.loc,'class')
                sortmatrix=horzcat(sortmatrix,obj.locData.loc.class);
            end
            sortmatrix=horzcat(sortmatrix,gval);
            [sm,indsort]=sortrows(sortmatrix); 
            listm=sm(:,end);
            list2=0*listm;
            ind=1;
            list2(1)=ind;
            for k=2:length(listm)
                if listm(k)~=listm(k-1)
                    ind=ind+1;
                end
                list2(k)=ind;
            end
            numbers=1:size(sortmatrix,1);
            indold=numbers(indsort);
            clear numbers
            [~,indback]=sort(indold);
            listback=list2(indback);
            %later:sortmatrix to not sort across files
             gr=Grouper(obj.locData);
             obj.locData.loc.groupindex=listback;
             gr.indsortlist=indsort;
%                
            gr.dllistsort=[diff(list2);1];
            gr.combine;
            
            obj.locData.filter;
             
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function pard=guidef(obj)
pard.locfield.object=struct('String','','Style','popupmenu');
pard.locfield.position=[3,2];
pard.locfield.Width=1;
pard.syncParameters={{'locFields','locfield',{'String'}}};
pard.plugininfo.description='Combines localizations which have the same value in a selected field';
pard.plugininfo.type='ProcessorPlugin';
end
