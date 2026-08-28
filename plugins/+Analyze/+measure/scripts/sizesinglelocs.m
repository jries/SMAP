% determine size of grouped localizations for filtering

% field ="tid"; % for MINFLUX
field="groupindex";
skipfirst=0; % for MINFLUX: skip to measure size without tails

valnan=100; % if only one localization is found, std would be zero, so set to a high value. Optional can also be set to zero to keep all single locs.

groupid=g.locData.loc.(field);
[sgroupid,indsort]=sort(groupid);
% groupids=unique(groupid);

xs=g.locData.loc.xnm(indsort);
ys=g.locData.loc.ynm(indsort);

stdxy=0*groupid;
ind1=1;
while ind1<length(sgroupid)
    ghere=sgroupid(ind1);
    ind2=ind1+1;
    while ind2<length(sgroupid) && sgroupid(ind2)==ghere
        ind2=ind2+1;
    end
    indh=ind1:ind2-1;
    if length(indh)<=skipfirst+1 
        stdxy(indh)=valnan;
    else
        stdxy(indh)=sqrt(std(xs(indh(skipfirst+1:end))).^2+std(ys(indh(skipfirst+1:end))).^2);
        % stdxy(indh)=mean(sgroupid(indh(skipfirst+1:end))); %test
    end
    ind1=ind2;
end

[~,ib]=sort(indsort);
stdxyb=stdxy(ib);
g.locData.loc.groupstd=stdxyb;
g.locData.regroup;
