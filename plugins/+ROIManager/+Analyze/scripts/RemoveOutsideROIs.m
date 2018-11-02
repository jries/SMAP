global se
rs=se.getPar('se_siteroi');
sites=se.sites;
locD=se.locData;
allind=false(length(locD.loc.xnm),1);
for k=1:length(sites)
    poss=sites(k).pos;
    loc=locD.getloc({'inungrouped'},'Position',[poss(1), poss(2), rs, rs],'filenumber',sites(k).info.filenumber);
    indh=loc.inungrouped;
    allind=allind | indh;
end
se.locData.removelocs(~allind);
se.locData.regroup;