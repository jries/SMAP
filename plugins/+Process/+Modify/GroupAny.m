classdef GroupAny<interfaces.DialogProcessor
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
%              gu=unique(gval);
%              mapping=zeros(max(gval),1);mapping(gu)=1:length(gu);
%              gvalm=mapping(gval);
             obj.locData.loc.groupindex=listback;
             gr.indsortlist=indsort;
%      
%              [sortedm,gr.indsortlist]=
%              listm=sortedm(:,end);
%              listsort
%             [listsort,gr.indsortlist]=sort(gval);
            
%             listm=mapping(listsort);
            
            %map to 
           
            
            
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
% pard.t1.object=struct('String','Only localizations displayed in layer 1 are used','Style','text');
% pard.t1.position=[1,1];
% pard.t1.Width=4;
% 
% pard.resultfield.object=struct('String','','Style','edit');
% pard.resultfield.position=[3,1];
% pard.resultfield.Width=0.8;
% 
% 
% pard.t2.object = struct('String','=','Style','text');
% pard.t2.position=[3,1.8];
% pard.t2.Width=0.2;

pard.locfield.object=struct('String','','Style','popupmenu');
pard.locfield.position=[3,2];
pard.locfield.Width=1;

% pard.filter.object=struct('String',{{'median','mean','max','min','std','5-95% percentile'}},'Style','popupmenu','Callback',{{@obj.locfield_callback}});
% pard.filter.position=[3,3];
% pard.filter.Width=1;
% 
% pard.t3.object = struct('String','Spatial scale (nm)','Style','text');
% pard.t3.position=[4,1];
% pard.t3.Width=1;
% 
% pard.spatialscale.object = struct('String','100','Style','edit');
% pard.spatialscale.position=[4,2];
% pard.spatialscale.Width=1;
% pard.t3.object=struct('String','Equation. Use fieldnames (e.g. xnm, phot) as variables','Style','text');
% pard.t3.position=[2,2];
% pard.t3.Width=3;
% 
% pard.equation.object = struct('String','(locprecnm<25 & PSFxnm>100) | numberInGroup>1','Style','edit');
% pard.equation.position=[3,2];
% pard.equation.Width=3;


% pard.dataselect.object=struct('Style','popupmenu','String','File');
% pard.dataselect.position=[2,1];
% pard.dataselect.object.TooltipString='choose localization file data set';
% 
% 
% pard.dataselectall.object=struct('Style','checkbox','String','all','Value',1);
% pard.dataselectall.position=[2,2];
% pard.dataselectall.object.TooltipString='apply on all';


pard.syncParameters={{'locFields','locfield',{'String'}}};

pard.plugininfo.type='ProcessorPlugin';
end
