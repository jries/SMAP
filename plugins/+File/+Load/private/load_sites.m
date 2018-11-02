function [locDataAll,sexp]=load_sites(filedat)


locDataAll=interfaces.LocalizationData;


filelist=filedat.siteparsave.filelist.list;
for k=1:length(filelist);
    f2.fitpos=filelist(k).pos.fitpos;
    f2.infox=filelist(k).par.info;
    f2.filename=filelist(k).par.file;
    [templocData,parameters,SEtemp]=recgui.loadfitposV2(f2);
    %put old filenumber in
    templocData.loc.filenumber= templocData.loc.filenumber*0+filelist(k).number;
    locDataAll.addLocData(templocData);
    templocData.files.file(k).number=filelist(k).number;
    locDataAll.files.file(k)=templocData.files.file(1);
    
end
locDataAll.files.filenumberEnd=max(locDataAll.loc.filenumber);


%sites
sexp=interfaces.SiteExplorer;
for k=1:length(filelist)
    file=filelist(k);
    sexp.addFile(file.filename,file.number,file.par.info)
end

celllist=filedat.siteparsave.celllist.list;
for k=1:length(celllist)   
    cell=celllist(k);
    cellnew=interfaces.SEsites;
    cellnew.pos=cell.pos*1000;
    cellnew.ID=cell.number;
    cellnew.info.filenumber=cell.file.number;
    sexp.cells(k)=cellnew;   
end

sitelist=filedat.siteparsave.sitelist.list;
for k=1:length(sitelist)   
    site=sitelist(k);
    sitenew=interfaces.SEsites;
    sitenew.pos=site.pos*1000;
    sitenew.ID=site.number;
    sitenew.info.filenumber=site.file.number;
    sitenew.info.cell=site.cell.number;
    sitenew.annotation=getannotation(site);   
    sexp.sites(k)=sitenew;   
end

function annotation=getannotation(site)
    %annotations
    siteinfo=site.siteinfo;
    for m=1:4
        annotation.(['list' num2str(m)]).value=siteinfo.parlist.values{m};
        annotation.(['list' num2str(m)]).string=siteinfo.parlist.text{m};
    end    
    annotation.comments=siteinfo.comments;   
    annotation.rotationpos=linefrompos(siteinfo.measure.linepos{3});
    annotation.line1=linefrompos(siteinfo.measure.linepos{1});
    annotation.line2=linefrompos(siteinfo.measure.linepos{2});

function line=linefrompos(pos)
pos2=pos;pos2(1,:)=pos(2,:);pos2(2,:)=pos(1,:);
line.pos=pos2;
line.length=sqrt(sum((pos2(1,:)-pos2(2,:)).^2))*1000;
line.angle=pos2angle(pos2);
line.value=line.length;

