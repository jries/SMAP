function bead=sites2beads(obj,p)

se=obj.locData.SE;
sites=se.sites;
for k=length(sites):-1:1
    bead(k).loc=obj.locData.getloc({'xnm','ynm','frame','PSFxnm','PSFynm','phot','znm'},'grouping','ungrouped','position',sites(k));
    bead(k).filenumber=sites(k).info.filenumber;
    bead(k).pos=sites(k).pos;   
    fn=fieldnames(bead(k).loc);
    for f=1:length(fn)
        bead(k).loc.(fn{f})=double(bead(k).loc.(fn{f}));
    end
end
use=logical(getFieldAsVector(sites,'annotation','use'));
bead=bead(use);
end