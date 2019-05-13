%import ROIs from Vilma and write attributes to localizations
addlocinfo=true;
addsite=false;

filename=g.locData.files.file(1).name;
[path,file,ext]=fileparts(filename);
infofile=[path filesep file '_3DVolume.txt'];
siteposfile=[path filesep file '_3DVolume_1.txt'];

classfile=[path filesep 'cluster_memberships-NUP-SNAP.txt'];
classtable=readtable(classfile);
for k=length(classtable.pore_source):-1:1    
    [~,ff]=fileparts(classtable.pore_source{k})  ;
    allclassfiles{k}=ff;
end

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
g.locData.loc.vilmaclass=g.locData.loc.xnm*0;
for k=1:length(siteattributes.id)
    pos=[siteattributes.y(k)*pixxy+posoff(1) siteattributes.x(k)*pixxy+posoff(2) siteattributes.z(k)*pixz+zpos];
    indc=siteattributes.id(k)==classtable.pore_id & contains(allclassfiles,file)';
    class=classtable.cluster_id((indc));
    if addsite
        currentsite=interfaces.SEsites;
        currentsite.pos=pos;     
        currentsite.info.cell=SE.currentcell.ID;
        currentsite.info.filenumber=SE.currentfile.ID;
        currentsite.evaluation.Vilma.siteid=siteattributes.id(k);
        currentsite.evaluation.Vilma.class=class;
        currentsite.ID=SE.addSite(currentsite);
    end
    if addlocinfo
        inh=(g.locData.loc.xnm-pos(1)).^2+(g.locData.loc.ynm-pos(2)).^2<R^2;
        g.locData.loc.vilmaid(inh)=siteattributes.id(k);
        if ~isempty(class)
        g.locData.loc.vilmaclass(inh)=class;
        end
    end
end

g.locData.regroup;
