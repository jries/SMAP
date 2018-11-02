function loco=get2Clocintensities(loc,transform,file,p)

 p.datapart.selection='all';
loct=apply_transform_locs(loc,transform,file,p);
indref=transform.getRef(loc.xnm,loc.ynm);
indreff=find(indref);
indtarf=find(~indref);

locr=copystructReduce(loc,indref);
loctr=copystructReduce(loct,~indref);

[iA,iB,uiA,uiB]=matchlocsall(renamefields(locr),renamefields(loctr),0,0,100);

iA=indreff(iA);
iB=indtarf(iB);
uiA=indreff(uiA);
uiB=indtarf(uiB);

loco.intA1=zeros(size(loc.xnm),'single')-100;
loco.intB1=zeros(size(loc.xnm),'single')-100;
loco.intA1(iA)=loc.phot(iA);
loco.intB1(iA)=loc.phot(iB);
loco.intA1(iB)=loc.phot(iA);
loco.intB1(iB)=loc.phot(iB);


loco.intA1(uiA)=loc.phot(uiA);
loco.intB1(uiB)=loc.phot(uiB);

function loco=renamefields(loci)
loco.x=loci.xnm;
loco.y=loci.ynm;
loco.frame=loci.frame;

