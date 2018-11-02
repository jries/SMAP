classdef SaveSiteParameters<interfaces.DialogProcessor&interfaces.SEProcessor
    properties
    end
    methods
        function obj=SaveSiteParameters(varargin)        
            obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={};
        end
        function initGui(obj)

        end
        
        function out=run(obj,p)  
            out=[];
           sitesin=obj.SE.sites;
           switch p.saveselection.selection
               case 'all site info'
                   fieldc={'pos','ID','info','annotation','evaluation','name'};
               case 'evaluation'
                   fieldc={'evaluation','name'};
           end
           
           for k=length(sitesin):-1:1
               sites(k)=copyfields([],sitesin(k),fieldc);
           end
           
           outf=obj.getPar('lastSMLFile');
           outf=strrep(outf,'_sml.mat','_ROIeval.mat');
           [file,pfad]=uiputfile(outf);
           if file
               save([pfad,file],'sites','-v7.3')
           end
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end


function pard=guidef(obj)
pard.saveselection.object=struct('String',{{'all site info','evaluation'}},'Style','popupmenu');
pard.saveselection.position=[1,1];
pard.saveselection.Width=2;

pard.plugininfo.type='ROI_Analyze';
end