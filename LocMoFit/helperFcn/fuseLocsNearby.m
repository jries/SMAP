function newLocs = fuseLocsNearby(locs, dist)
    neighbors = rangesearch([locs.xnm,locs.ynm,locs.znm], [locs.xnm,locs.ynm,locs.znm],dist);
    allGroups = {};
    neighborsMatrix = false(length(neighbors));
    for k = 1:length(neighbors)
        idxNeighbors = neighbors{k};
        neighborsMatrix(k,idxNeighbors)=true;
    end
    neighborsMatrix = sparse(neighborsMatrix);
    nMember = sum(neighborsMatrix,2);
    idxRemain = 1:length(neighbors);
    idxRm = 1:length(neighbors);
    while any(nMember>1)
        [~,idxTakeOut] = max(nMember);
        idxMember = neighborsMatrix(idxTakeOut,:);
        allGroups{end+1} = idxRemain(idxMember);
        idxRemain(idxMember) = [];
        neighborsMatrix(idxMember,:)=[];
        neighborsMatrix(:,idxMember)=[];
        nMember = sum(neighborsMatrix,2);
    end
    locs.n = ones(size(locs.xnm));
    for k = 1:length(allGroups)
        oneMembers = allGroups{k};
        locs.xnm(end+1) = mean(locs.xnm(oneMembers));
        locs.ynm(end+1) = mean(locs.ynm(oneMembers));
        locs.znm(end+1) = mean(locs.znm(oneMembers));
        locs.n(end+1) = length(oneMembers);
    end
    lIdxRm = ~ismember(idxRm, idxRemain);
    idxRm = idxRm(lIdxRm);
    newLocs.xnm = locs.xnm;
    newLocs.ynm = locs.ynm;
    newLocs.znm = locs.znm;
    newLocs.n = locs.n;
    newLocs.layer = ones(size(locs.xnm));
    fn = fieldnames(newLocs);
    for k = 1:length(fn)
        newLocs.(fn{k})(idxRm) = [];
    end
end
