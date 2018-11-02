classdef ShowHistory<interfaces.DialogProcessor
    methods
        function obj=ShowHistory(varargin)    
            obj@interfaces.DialogProcessor(varargin{:}) ;
        end
        
        function out=run(obj,p)
            hist=obj.locData.history;
            texta={};
            
            ttc=[];
            for k=1:length(hist)
                
                phist=hist{k};
                txt=struct2txt(phist,'');
                if isfield(phist,'name')
                    name=phist.name;
                else 
                    name='';
                end
                texta{end+1}=['Module' num2str(k) ': ' name];
                texta(end+1:end+length(txt))=txt;
                texta{end+1}='';
                
                
                ttc=[ttc 'Module' num2str(k) ': ' name 10 ];
                for l=1:length(txt)
                    
                    th=strrep(txt{l},'=',['=' 9]);
                    ttc=[ttc th 10];
                end
                
            end
            
            listdlg('ListString',texta,'ListSize',[800,800]);
            out=texta;
            
      
            clipboard('copy',ttc)
           
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        

    end
end




function pard=guidef
pard.plugininfo.name='Show History';
pard.plugininfo.type='ProcessorPlugin';
end