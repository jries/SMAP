function locs = getLocsFromSite(g, siteID, rmFilter, center2Ori)
    se = g.locData.SE;
    sites = g.locData.SE.sites;
    roiSize = g.getPar('se_siteroi');
    ID = getFieldAsVector(sites, 'ID');
    [~,idxSite] = ismember(siteID, ID);
    k = 1;
    % hack an evaluate plug-in in order to use the obj.getLocs(...)
    fdcal=figure('Name','temp', 'Visible', 'off');
    dcal=plugin('ROIManager','Evaluate','generalStatistics',fdcal,g.P);
    dcal.attachLocData(se.locData);
    dcal.makeGui;
    dcal.site=sites(idxSite(k));
    dcal.site.image = se.plotsite(sites(k));
    [locs,indlocOne] = dcal.getLocs({'xnmrot','ynmrot','znm','locprecnm', 'locprecznm','channel'},'size',roiSize,'grouping', 'grouped'); % per ROI info.
    fn = fieldnames(locs);
    lEmpty = structfun(@isempty, locs);
    locs = rmfield(locs, fn(lEmpty));
    close(fdcal)
end