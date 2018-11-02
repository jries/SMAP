function bead=segmentb(obj,p)
if isfield(p,'zrangeuse')
minframes=(p.zrangeuse(2)-p.zrangeuse(1))/ p.dz/2;
else
    minframes=800/p.dz;
end
% minframes=(p.zrangeuse(2)-p.zrangeuse(1))/ p.dz/5;

locDatacopy=obj.locData.copy;
locDatacopy.regroup(150,minframes/2);


% beadind=find(locDatacopy.grouploc.numberInGroup>minframes);
% xb=locDatacopy.grouploc.xnm(beadind);
% yb=locDatacopy.grouploc.ynm(beadind);
% filenumber=locDatacopy.grouploc.filenumber(beadind);
locg=locDatacopy.getloc({'xnm','ynm','numberInGroup','filenumber'},'layer',1,'Position','roi','removeFilter','filenumber','grouping','grouped');
beadind=find(locg.numberInGroup>minframes);
xb=locg.xnm(beadind);
yb=locg.ynm(beadind);
filenumber=locg.filenumber(beadind);
winsize=250;
for k=length(beadind):-1:1
%     beadhere=(mywithin(locDatacopy.loc.xnm,[xb(k)-winsize/2 winsize],locDatacopy.loc.ynm,[yb(k)-winsize/2 winsize]));
%     beadfile=locDatacopy.loc.filenumber==filenumber(k);
%     beadhere=beadhere&beadfile;
    bead(k).loc=locDatacopy.getloc({'xnm','ynm','PSFxnm','PSFynm','frame','phot','znm'},'filenumber',filenumber(k),'Position',[xb(k) yb(k) winsize]);
%     bead(k).loc=copystructReduce(locDatacopy.loc,beadhere,{'xnm','ynm','PSFxnm,','PSFynm','frame'});
    bead(k).filenumber=filenumber(k);
    bead(k).pos=[xb(k) yb(k)];
    fn=fieldnames(bead(k).loc);
    for f=1:length(fn)
        bead(k).loc.(fn{f})=double(bead(k).loc.(fn{f}));
    end    
end

%remove beads with close neighbours.
mindist=1000;
goodind=true(length(beadind),1);
for k=1:length(beadind)
    for l=k+1:length(beadind)
        d2=sum((bead(k).pos-bead(l).pos).^2);
        df=mean(bead(k).loc.frame)-mean(bead(l).loc.frame);
        if bead(k).filenumber==bead(l).filenumber && d2<mindist^2 %&& df<p.zrangeuse(2)-p.zrangeuse(1)
            goodind(k)=false;
            goodind(l)=false;
        end
    end
end
bead=bead(goodind);
end
