function newLocs = fuseLocsNearby(locs, dist)
useNew = 1;
if useNew
    neighbors = rangesearch([locs.xnm,locs.ynm,locs.znm], [locs.xnm,locs.ynm,locs.znm],dist);
    numOfNB = cellfun(@length, neighbors);
    
    connections = zeros(sum(numOfNB-1)/2,2);
    currentS = 0;
    for k = 1:length(neighbors)
        idx = neighbors{k}>k;
        nn = sum(idx);
        if nn>0
            connections(currentS+1:currentS+nn,1) = k;
            connections(currentS+1:currentS+nn,2) = neighbors{k}(idx);
        end
        currentS = currentS+nn;
    end
    
    grpIdx = zeros(length(neighbors),1);
    currentGrp = 0;
    while any(numOfNB>1)
        currentGrp = currentGrp+1;
        [~,idxTakeOut] = max(numOfNB);
        selectedCon = connections(connections(:,1)==idxTakeOut|connections(:,2)==idxTakeOut,:);
        primaryNB = unique(selectedCon(:));
        grpIdx(primaryNB) = currentGrp;
        con2rm = ismember(connections,primaryNB);
        lCon2rm = sum(con2rm,2)>=1;
        con2rm = connections(lCon2rm,:);
        nNBRm = histcounts(con2rm(:),1:length(numOfNB)+1);
        numOfNB = numOfNB - nNBRm';
        connections(lCon2rm,:)=[];
    end
    currentGrp = currentGrp+1;
    lZero = grpIdx==0;
    grpIdx(lZero) = currentGrp:currentGrp-1+sum(lZero);
    nLocsFused = accumarray(grpIdx, 1);
    newLocs.xnm = accumarray(grpIdx, locs.xnm)./nLocsFused;
    newLocs.ynm = accumarray(grpIdx, locs.ynm)./nLocsFused;
    newLocs.znm = accumarray(grpIdx, locs.znm)./nLocsFused;
    newLocs.n = nLocsFused;
    newLocs.layer = ones(size(newLocs.xnm));
else
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
end
