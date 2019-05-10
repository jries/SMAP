%import ROIs from Vilma and write attributes to localizations
filename=g.locData.files.file(1).name;
[path,file,ext]=fileparts(filename);
infofile=[path filesep file '_3DVolume.txt'];
siteposfile=[path filesep file '_3DVolume_1.txt'];

siteattributes=readtable(siteposfile);
infot=readtable(infofile,'ReadVariableNames',false,'ReadRowNames',true);
sr_pos=str2num(infot.Var1{6});
sr_size=str2num(infot.Var1{7});
pixxy=str2num(infot.Var1{8});
pixz=str2num(infot.Var1{9});
minz=str2num(infot.Var1{4});
maxz=str2num(infot.Var1{5});
zpos=mean(minz,maxz);

posoff=sr_pos-sr_size;
SE=g.locData.SE;
R=g.getPar('se_siteroi')/2;
g.locData.loc.vilmaid=g.locData.loc.xnm*0;
for k=1:length(siteattributes.id)
    pos=[siteattributes.y(k)*pixxy+posoff(1) siteattributes.x(k)*pixxy+posoff(2) siteattributes.z(k)*pixz+zpos];
    currentsite=interfaces.SEsites;
    currentsite.pos=pos;     
    currentsite.info.cell=SE.currentcell.ID;
    currentsite.info.filenumber=SE.currentfile.ID;
    currentsite.evaluation.Vilma.siteid=siteattributes.id(k);
    currentsite.ID=SE.addSite(currentsite);
    inh=(g.locData.loc.xnm-pos(1)).^2+(g.locData.loc.ynm-pos(2)).^2<R^2;
    g.locData.loc.vilmaid(inh)=siteattributes.id(k);
end

g.locData.regroup;
