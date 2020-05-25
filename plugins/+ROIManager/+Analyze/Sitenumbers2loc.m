classdef Sitenumbers2loc<interfaces.DialogProcessor&interfaces.SEProcessor
%     Adds two field to the localization data containing the site number
%     and the cell number, respectively.
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
            roisize=obj.getPar('se_siteroi');
            fovsize=obj.getPar('se_sitefov');
            for k=1:length(sites)
                posxy=sites(k).pos(1:2);
                switch p.roimode.selection   
                    case 'site FoV'
                        position=[posxy fovsize fovsize];
                    case 'ROI square'
                        position=[posxy roisize roisize];
                    case 'ROI round'
                        position=[posxy roisize/2];
                end
                [l,ind]=obj.locData.getloc({'filenumber','xnm','ynm'},'Position',position);
%                 find=find(ind);
%                 indfile=l.filenumber==sites(k).info.filenumber;
%                 indh=find(indfile);
% figure(88);plot(l.xnm,l.ynm,'.')
                indh=ind & (obj.locData.loc.filenumber==sites(k).info.filenumber);
                sitenumbers(indh)=sites(k).ID;
                cellnumbers(indh)=sites(k).info.cell;
            end
            obj.locData.setloc('sitenumbers',sitenumbers);
            obj.locData.setloc('cellnumbers',cellnumbers);
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

pard.roimode.object=struct('String',{{'site FoV','ROI square','ROI round'}},'Style','popupmenu');
pard.roimode.position=[2,1];
pard.roimode.Width=1;


pard.plugininfo.description='Adds two field to the localization data containing the site number and the cell number, respectively.';
pard.plugininfo.type='ROI_Analyze';
end