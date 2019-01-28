classdef Sitenumbers2loc<interfaces.DialogProcessor&interfaces.SEProcessor
    methods
        function obj=Sitenumbers2loc(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
%             obj.inputParameters={'se_layerson'};
            obj.showresults=false;
        end
        
        function out=run(obj,p)  
            out=[];
            obj.setPar('undoModule','Sitenumbers2loc');
            notify(obj.P,'backup4undo');
            
            sites=obj.SE.sites;
            
            ld=obj.locData.loc;
            sitenumbers=0*ld.xnm;
            cellnumbers=0*ld.xnm;
            for k=1:length(sites)
                [l,ind]=obj.locData.getloc('filenumber','Position',sites(k));
%                 find=find(ind);
%                 indfile=l.filenumber==sites(k).info.filenumber;
%                 indh=find(indfile);
                indh=ind & (obj.locData.loc.filenumber==sites(k).info.filenumber);
                sitenumbers(indh)=sites(k).ID;
                cellnumbers(indh)=sites(k).info.cell;
            end
            obj.locData.addloc('sitenumbers',sitenumbers);
            obj.locData.addloc('cellnumbers',cellnumbers);
            obj.locData.regroup;
          
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end




function pard=guidef

pard.t1.object=struct('String','Write site and cell numbers to locData.loc','Style','text');
pard.t1.position=[1,1];
pard.t1.Width=4;


pard.plugininfo.type='ROI_Analyze';


end