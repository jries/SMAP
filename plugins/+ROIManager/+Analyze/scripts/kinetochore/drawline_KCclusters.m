%generate line along spindels for double selected kinetochore clusters
%through the center of the selected ROIs
% cellline='line2';
linenumber=2;

sites=g.locData.SE.sites;
cells=getFieldAsVector(sites,'info.cell');
cellnumbers=unique(cells);
% annotatedsites=false(length(sites),1);

for k=1:length(cellnumbers)
    incell=find(cells==cellnumbers(k));
    if length(incell)~=2
        disp(['cell ' num2str(cellnumbers(k)) ' does not have two sites'])
        continue
    end
    pos1=sites(incell(1)).pos(1:2);
    pos2=sites(incell(2)).pos(1:2);
    
    dpos1=(pos2-pos1);
    len1=sqrt(sum(dpos1.^2));
    dpos1n=dpos1/len1*g.getPar('se_siteroi');
    linepos1=[-dpos1n/2; dpos1n/2]+sites(incell(1)).pos(1:2);
    
%     sites(incell(1)).annotation.(celline).pos=linepos1/1000;
    sites(incell(1)).setlinepos(linenumber,linepos1)
    
    linepos1=[dpos1n/2; -dpos1n/2]+sites(incell(2)).pos(1:2);
    sites(incell(2)).setlinepos(linenumber,linepos1)
    
%     if linestruct.length>0
%         samecell=find(cells==sites(k).info.cell);
%         annotatedsites(k)=true;
%         for s=1:length(samecell)
%             if ~annotatedsites(samecell(s))
%                 pos=linestruct.pos;
%                 newpos=sites(samecell(s)).pos(1:2)/1000;
%                 pos=pos-mean(pos,1)+newpos;
%                 sites(samecell(s)).annotation.(celline)=linestruct;
%                 sites(samecell(s)).annotation.(celline).pos=pos;
%             end
%         end
%     end
end
            