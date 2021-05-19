% plot drift of beads
locs=g.locData.getloc({'xnm','ynm','numberInGroup','groupindex','frame'},'Position','roi','grouping','ungrouped');
groups=unique(locs.groupindex);
nmax=max(locs.numberInGroup);

indgood=true(length(locs.xnm),1);
for k=1:length(groups)
   ingroup=locs.groupindex == groups(k);
   if mean(locs.numberInGroup(ingroup)) < nmax
       indgood(ingroup)=false;
   end
end
nframes=max(locs.frame);
xm=zeros(nframes,1);
ym=zeros(nframes,1);
for k=1:nframes
    inframe=locs.frame==k;
    inh=inframe&indgood;
    xm(k)=mean(locs.xnm(inh));
    ym(k)=mean(locs.ynm(inh));
end

frames=(1:nframes)';
figure(88);
hold on
plot(frames,xm-min(xm),frames,ym-min(ym))
xlabel('frame')
ylabel('x / y (nm)');
legend('x','y')