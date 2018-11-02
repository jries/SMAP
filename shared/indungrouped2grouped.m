function ind=indungrouped2grouped(indcombined,groupindex)
gind=unique(groupindex(indcombined));
ind=false(length(indcombined),1);
ind(gind)=true;
end