function bead=segmentb_so(locs,dz)
minframes=800/dz;
sortmatrix=horzcat(locs.filenumber,locs.frame,locs.x);
[~,indsort]=sortrows(sortmatrix,1:3);
list=connectsingle2c(double(locs.x(indsort)),double(locs.y(indsort)),double(locs.frame(indsort)),double(0.5),int32(0),int32(10000));
hh=histcounts(list,1:max(list)+1);
beadind=find(hh>=minframes);
winsize=2.5;
for k=length(beadind):-1:1
    indh=list==beadind(k);indh1=find(indh,1,'first');
    beadposx=median(locs.x(indh));
    beadposy=median(locs.y(indh));
    beadfile=locs.filenumber(indh1);
    inrange=(locs.x-beadposx).^2+(locs.y-beadposy).^2<winsize^2&locs.filenumber==beadfile;
    bead(k).filenumber=beadfile;
    bead(k).pos=[beadposx beadposy];
    bead(k).loc.z=double(locs.z(inrange));
    bead(k).loc.phot=double(locs.phot(inrange));
    bead(k).loc.frame=double(locs.frame(inrange)); 
end

%remove beads with close neighbours.
mindist=10;
goodind=true(length(beadind),1);
for k=1:length(beadind)
    for l=k+1:length(beadind)
        d2=sum((bead(k).pos-bead(l).pos).^2);
        if bead(k).filenumber==bead(l).filenumber && d2<mindist^2
            goodind(k)=false;
            goodind(l)=false;
        end
    end
end
bead=bead(goodind);
end
