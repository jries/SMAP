classdef FilterSites<interfaces.DialogProcessor&interfaces.SEProcessor
    properties
    end
    methods
        function obj=FilterSites(varargin)        
            obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={};
        end
        function initGui(obj)

        end
        
        function out=run(obj,p)  
            out=[];
            for e=1:3
                eq{e}=p.(['eq' num2str(e)]);
                for k=1:6
                    if isfield(p,['x' num2str(k)])
                        eq{e}=strrep(eq{e},['x' num2str(k)],['sites(siteind).' p.(['x' num2str(k)])]);
                    end
                end     
            end
            sites=obj.SE.sites;
            eqtot=['(' eq{1}];
            
            for e=2:3
                if strcmp(p.(['logiceq' num2str(e)]).selection,'AND')
                    operator='&';
                else
                    operator='|';
                end
                if any(strfind(eq{e},'sites'))
                
                eqtot=[eqtot ')' operator '(' eq{e}];
                end
            end
            eqtot=[eqtot ')'];
            disp(eqtot)
            for siteind=1:length(sites)
                sites(siteind).annotation.use=eval(eqtot);
            end
            obj.SE.processors.preview.updateSitelist;
   
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function xselect(obj,callobj,data,field)
            site=obj.SE.currentsite;
            siteh=site;
            str=obj.guihandles.(field).String;
            ind=strfind(str,'.');
            oldstr='';
            if length(ind)>1
                oldstr=str(1:ind(end-1)-1);
                siteh=eval(['site.' oldstr]);
                oldstr=[oldstr '.'];
                str=str(ind(end-1)+1:ind(end)-1);
            elseif length(ind)==1
                str=str(1:ind(end)-1);
            end
            
%             val=eval(['site.' str]);
%             if ~isstruct(val)
%                 ind=strfind(str,'.');
%                 str=str(1:ind(end)-1)
%                 val=eval(['site.' str]);
%             end
            stro=browsefields(siteh,str);
            obj.guihandles.(field).String=[oldstr stro];

        end
        function useall(obj,a,b)
            sites=obj.SE.sites;
            for siteind=1:length(sites)
                
                sites(siteind).annotation.use=true;
                obj.SE.processors.preview.updateSitelist;
            end
        end
    end
end


function pard=guidef(obj)

pard.x1b.object=struct('String','x1','Style','pushbutton','Callback',{{@obj.xselect,'x1'}});
pard.x1b.position=[1,1];
pard.x1b.Width=0.2;

pard.x1.object=struct('String','evaluation','Style','edit');
pard.x1.position=[1,1.2];
pard.x1.Width=1.8;

pard.x2b.object=struct('String','x2','Style','pushbutton','Callback',{{@obj.xselect,'x2'}});
pard.x2b.position=[1,3];
pard.x2b.Width=0.2;

pard.x2.object=struct('String','evaluation','Style','edit');
pard.x2.position=[1,3.2];
pard.x2.Width=1.8;

pard.x3b.object=struct('String','x3','Style','pushbutton','Callback',{{@obj.xselect,'x3'}});
pard.x3b.position=[2,1];
pard.x3b.Width=0.2;

pard.x3.object=struct('String','','Style','edit');
pard.x3.position=[2,1.2];
pard.x3.Width=1.8;

pard.x4b.object=struct('String','x4','Style','pushbutton','Callback',{{@obj.xselect,'x4'}});
pard.x4b.position=[2,3];
pard.x4b.Width=0.2;

pard.x4.object=struct('String','','Style','edit');
pard.x4.position=[2,3.2];
pard.x4.Width=1.8;

pard.x5b.object=struct('String','x5','Style','pushbutton','Callback',{{@obj.xselect,'x5'}});
pard.x5b.position=[3,1];
pard.x5b.Width=0.2;

pard.x5.object=struct('String','annotation','Style','edit');
pard.x5.position=[3,1.2];
pard.x5.Width=1.8;

pard.x6b.object=struct('String','x6','Style','pushbutton','Callback',{{@obj.xselect,'x6'}});
pard.x6b.position=[3,3];
pard.x6b.Width=0.2;

pard.x6.object=struct('String','annotation','Style','edit');
pard.x6.position=[3,3.2];
pard.x6.Width=1.8;

pard.eq1.object=struct('String',' ','Style','edit');
pard.eq1.position=[4,1.7];
pard.eq1.Width=3.3;

pard.eq2.object=struct('String',' ','Style','edit');
pard.eq2.position=[5,1.7];
pard.eq2.Width=3.3;

pard.logiceq2.object=struct('String',{{'AND','OR'}},'Style','popupmenu');
pard.logiceq2.position=[5,1];
pard.logiceq2.Width=.7;

pard.eq3.object=struct('String',' ','Style','edit');
pard.eq3.position=[6,1.7];
pard.eq3.Width=3.3;

pard.logiceq3.object=struct('String',{{'AND','OR'}},'Style','popupmenu');
pard.logiceq3.position=[6,1];
pard.logiceq3.Width=.7;

pard.reset.object=struct('String','use all','Style','pushbutton','Callback',{{@obj.useall}});
pard.reset.position=[8,1];
pard.reset.Width=1;

pard.plugininfo.type='ROI_Analyze';
end