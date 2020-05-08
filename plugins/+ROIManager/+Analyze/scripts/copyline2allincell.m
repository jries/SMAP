%script to copy line to all other sites in same cell
celline='line2';

sites=g.locData.SE.sites;
cells=getFieldAsVector(sites,'info.cell');
annotatedsites=false(length(sites),1);

for k=1:length(sites)
    linestruct=sites(k).annotation.(celline);
    if linestruct.length>0
        samecell=find(cells==sites(k).info.cell);
        annotatedsites(k)=true;
        for s=1:length(samecell)
            if ~annotatedsites(samecell(s))
                pos=linestruct.pos;
                newpos=sites(samecell(s)).pos(1:2)/1000;
                pos=pos-mean(pos,1)+newpos;
                sites(samecell(s)).annotation.(celline)=linestruct;
                sites(samecell(s)).annotation.(celline).pos=pos;
            end
        end
    end
end
            