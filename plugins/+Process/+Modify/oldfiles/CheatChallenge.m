classdef CheatChallenge<interfaces.DialogProcessor
    methods
        function obj=CheatChallenge(varargin)   
            obj@interfaces.DialogProcessor(varargin{:});  
        end
        
        function out=run(obj,p)
            out=[];
            obj.setPar('undoModule','CheatChallenge');
            notify(obj.P,'backup4undo');
            
            locs=obj.locData.getloc({'xnm','ynm','znm','locprecnm','locprecznm'},'position','all');
            if isempty(locs.znm)
                locs.znm=0*locs.xnm;
                locs.locprecznm=locs.locprecnm;
            end
            xrange=[min(locs.xnm),max(locs.xnm)];
            yrange=[min(locs.ynm),max(locs.ynm)];
%             img=myhist2(locs.xnm,locs.ynm,px,px,xrange,yrange);
            [xs,sortind]=sort(locs.xnm);
            ys=locs.ynm(sortind);
%             if ~isempty(locs.znm)
            zs=locs.znm(sortind);
%             end
            ind1=1;
%             dxy=40;
%             dz=50;
            lx=length(xs);
            rmt=p.mtradius;
            lenfac=p.searchscale;
            minneighbours=p.minneigh;
            for k=1:lx
                sxy=locs.locprecnm(k);
                sz=locs.locprecznm(k);
                dxy=(1*2*sxy+rmt*1)*lenfac;
                dz=(1*sqrt(2)*sz+rmt*1)*lenfac;
                xh=xs(k);
                while(xs(ind1)<xh-dxy)&& ind1<lx ind1=ind1+1; end
                ind2=ind1;
                while(xs(ind2)<xh+dxy) && ind2<lx
                    ind2=ind2+1; 
                end

                distx=xs(ind1:ind2)-xs(k);
                disty=ys(ind1:ind2)-ys(k);
                distz=zs(ind1:ind2)-zs(k);
                indin=find((disty).^2/dxy^2 + (distx).^2/dxy^2+(distz).^2/dz^2<1);
                if length(indin)<minneighbours
                    continue
                end
                mx=mean(distx(indin));my=mean(disty(indin));mz=mean(distz(indin));
%                 mx=median(distx(indin));my=median(disty(indin));mz=median(distz(indin));
                %correct xy
                
                maxlocp=sxy*3;
                
                dvec=[mx, my, mz];
                
                dvecmt=dvec/norm(dvec)*rmt;
                dvecc=(dvec-1*dvecmt)*1;
                if norm(dvecc)>sqrt(2)*sxy;
                    dvecc=dvecc/norm(dvecc)*sqrt(2)*sxy;
                end
%                 maximum shift: sqrt(2)*locprec
                
                
                if  (mx^2+my^2)>rmt^2 && mx^2+my^2<(rmt+maxlocp)^2
                    locs.xnm(sortind(k))=locs.xnm(sortind(k))+dvecc(1);
                    locs.ynm(sortind(k))=locs.ynm(sortind(k))+dvecc(2);         
                end
                

                
                
%                 correct z
                maxlocpz=sz*4;
                if mz^2>rmt^2 && mz^2<(rmt+maxlocpz)^2
                    loc.znm(sortind(k))=loc.znm(sortind(k))+dvecc(3);
%                     loc.ynm(sortind(k))=loc.ynm(sortind(k))+my;         
                end
                %somewhere: check for number of neighbours: don't correct
                %background noise.
                
                %later: calculate density and correction with grouped data.
%                 dist=(disty).^2/dxy^2 + (distx).^2/dxy^2+(zs(ind1:ind2)-zs(k))^2/dz^2;
%                 dist=(ys(ind1:ind2)-ys(k)).^2/dxy^2 + (xs(ind1:ind2)-xs(k)).^2/dxy^2+(distz)^2/dz^2;
%                 indin=find(dist<1);
%                 indin=indin+ind1-1;
                
            end
                            obj.locData.loc.xnm=locs.xnm;
                obj.locData.loc.ynm=locs.ynm;
            obj.locData.regroup;   
        end
        function pard=guidef(obj)
            pard=guidef;
        end     

    end
end




function pard=guidef
pard.textb.object=struct('String','Sharpen Microtubul features','Style','text');
pard.textb.position=[1,1];
pard.textb.Width=3;
pard.textb.object.TooltipString='';

pard.mtradiust.object=struct('String','Radius MT (nm)','Style','text');
pard.mtradiust.position=[2,1];

pard.mtradius.object=struct('String','25','Style','edit');
pard.mtradius.position=[2,2];

pard.minneight.object=struct('String','min # neighbours','Style','text');
pard.minneight.position=[3,1];

pard.minneigh.object=struct('String','10','Style','edit');
pard.minneigh.position=[3,2];

pard.searchscalet.object=struct('String','factor ROI','Style','text');
pard.searchscalet.position=[4,1];

pard.searchscale.object=struct('String','1','Style','edit');
pard.searchscale.position=[4,2];
pard.plugininfo.type='ProcessorPlugin';
end