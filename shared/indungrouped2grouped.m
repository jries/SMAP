function ind=indungrouped2grouped(indcombined,groupindex)
gind=unique(groupindex(indcombined),'stable');
ind=false(max(groupindex),1);
ind(gind)=true;
end