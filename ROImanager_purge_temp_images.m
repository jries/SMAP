%purge images
se=g.locData.SE;

for k=1:length(se.sites)
    se.sites(k).image=[];
end
for k=1:length(se.cells)
    se.cells(k).image=[];
end
for k=1:length(se.files)
    se.files(k).image=[];
end