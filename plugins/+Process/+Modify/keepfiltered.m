classdef keepfiltered<interfaces.DialogProcessor
%     Keeps only filtered localizations. Can be used to reduce file size.
    methods
        function obj=keepfiltered(varargin)     
            obj@interfaces.DialogProcessor(varargin{:}) ;  
        end
        
        function out=run(obj,p) 
            out=[];
                obj.setPar('undoModule','RemoveLocs');
                notify(obj.P,'backup4undo');
                
                layers=find(obj.getPar('sr_layerson'));
                [l,indg]=obj.locData.getloc({'inungrouped'},'layer',layers(1),'Position','all','removeFilter',{'filenumber'});
                indgood=l.inungrouped;
                for k=2:length(layers)
                     [l,indg]=obj.locData.getloc({'inungrouped'},'layer',layers(k),'Position','all','removeFilter',{'filenumber'});
                     indgood=indgood | l.inungrouped;
                end
                indbad=~indgood;
                obj.locData.removelocs(indbad)
                obj.locData.filter;
                obj.locData.regroup; 
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end




function pard=guidef



pard.textb.object=struct('String','keeps only the filtered localizations','Style','text');
pard.textb.position=[1,1];
pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description='Keeps only filtered localizations. Can be used to reduce file size.';
end