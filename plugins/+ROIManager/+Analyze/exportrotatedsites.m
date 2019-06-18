sites=g.locData.SE.sites;
dy=300; %width of profile

dind=10000;
out=zeros(dind,7);
ind=1;
for k=1:length(sites)

layers=find(g.getPar('sr_layerson'));
locs=g.locData.getloc({'xnm','ynm','znm','channel','frame','locprecnm'},'layer',layers,'Position',sites(k));
pos=sites(k).annotation.rotationpos.pos;
angle=sites(k).annotation.rotationpos.angle;
len=sites(k).annotation.rotationpos.length;
if len==0
    continue
end
xy=mean(pos,1)*1000;
xr=cosd(angle)*(locs.xnm-xy(1))+sind(angle)*(locs.ynm-xy(2));
yr=cosd(angle)*(locs.ynm-xy(2))-sind(angle)*(locs.xnm-xy(1));

indg=abs(xr)<len & abs(yr) < dy/2;

numlocs=sum(indg);
if size(out,2)<ind+numlocs
    out(ind+dind,1)=0;
end
out(ind:ind+numlocs-1,:)=horzcat(xr(indg),yr(indg),locs.znm(indg),locs.channel(indg),...
    locs.frame(indg),locs.locprecnm(indg), k+0*xr(indg));
ind=ind+numlocs;
end
out(ind:end,:)=[];
to=array2table(out,'VariableNames',{'xnm','ynm','znm','channel','frame','locprecnm','ID'});

outf=[g.locData.files.file(1).name(1:end-8) '_export.mat'];
[f,path]=uiputfile(outf);
if f
    siteexport=out;
    siteexport=to;
    save([path f],'siteexport')
end