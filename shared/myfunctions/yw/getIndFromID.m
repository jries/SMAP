function idxSite = getIndFromID(g, siteID)
    sites = g.locData.SE.sites;
    ID = getFieldAsVector(sites, 'ID');
    [~,idxSite] = ismember(siteID, ID);
end