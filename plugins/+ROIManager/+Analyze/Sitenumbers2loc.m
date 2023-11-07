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
            locMoFitPar = sitenumbers;
            siteorder=0*ld.xnm;
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
                % look for overlapping sites, then use closest 
                assignedind=(sitenumbers>0) & indh;  %later populate manually if too slow
                
                if p.lLocMoFitPar
                    mp = 12;
                    parName = 'ringDistance';
                end

                if any(assignedind)
                    newindall=find(indh);
                    oldsites=sitenumbers(assignedind);
                    usid=unique(oldsites);
                    posnh=sites(k).pos;
                    for kh=1:length(usid)
                        idh=obj.locData.SE.indexFromID(obj.locData.SE.sites,usid(kh));
                        posoldh=sites(idh).pos;
                       
                        dold=(l.xnm-posoldh(1)).^2+(l.ynm-posoldh(2)).^2;
                        dnew=(l.xnm-posnh(1)).^2+(l.ynm-posnh(2)).^2;
                        indold=dold<dnew;
                        newindall(indold)=0;
                    end
                    newindall(newindall==0)=[];
                    sitenumbers(newindall)=sites(k).ID;
                    siteorder(newindall)=k;
                    cellnumbers(newindall)=sites(k).info.cell;
                    locMoFitPar(newindall)=sites(k).evaluation.LocMoFitGUI_3.allParsArg.value(mp);
                else
                    sitenumbers(indh)=sites(k).ID;
                    siteorder(indh)=k;
                    cellnumbers(indh)=sites(k).info.cell;
                    locMoFitPar(indh)=sites(k).evaluation.LocMoFitGUI_3.allParsArg.value(mp);
                end


                

            end
            if p.lSiteNumber
                obj.locData.setloc('sitenumbers',sitenumbers);
            end
            if p.lCellNumber
                obj.locData.setloc('cellnumbers',cellnumbers);
            end
            if p.lSiteOrder
                obj.locData.setloc('siteorder',siteorder);
            end
            if p.lLocMoFitPar
                obj.locData.setloc([parName '_LocMoFitGUI_3'],locMoFitPar);
            end
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

pard.roimode.object=struct('String',{{'site FoV','ROI square','ROI round'}},'Style','popupmenu','Value',3);
pard.roimode.position=[2,1];
pard.roimode.Width=1;

pard.lSiteNumber.object=struct('String','site number (ID)','value',1,'Style','checkbox');
pard.lSiteNumber.position=[3,1];
pard.lSiteNumber.Width=1;

pard.lSiteOrder.object=struct('String','site order','value',0, 'Style','checkbox');
pard.lSiteOrder.position=[3,2];
pard.lSiteOrder.Width=1;

pard.lCellNumber.object=struct('String','cell number','value',1,'Style','checkbox');
pard.lCellNumber.position=[3,3];
pard.lCellNumber.Width=1;

pard.lLocMoFitPar.object=struct('String','LocMoFit','value',0,'Style','checkbox');
pard.lLocMoFitPar.position=[3,4];
pard.lLocMoFitPar.Width=1;

pard.plugininfo.description='Adds two field to the localization data containing the site number and the cell number, respectively.';
pard.plugininfo.type='ROI_Analyze';
end