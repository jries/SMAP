function newLocs = groupLocs(locs, minDis)
    locs.uniqueID = (1:length(locs.xnm))';
    if ~isfield(locs, 'znm')
        pos = [locs.xnm locs.ynm];
        fn = {'xnm','ynm','layer'};
    else
        pos = [locs.xnm locs.ynm locs.znm];
        fn = {'xnm','ynm','znm','layer'};
    end
    allPairs = rangesearch(pos,pos,minDis);
    A = zeros(length(locs.xnm));
    for k = 1:length(locs.xnm)
        idxLinked = allPairs{k};
        if length(idxLinked)>1
            for l = 1:length(idxLinked)
                if k ~=idxLinked(l)
                    A(k,idxLinked(l))=1;
                end
            end
        end
    end
    nNeighbors = sum(A);
    AA = A(:,nNeighbors>0);
    AA = AA(nNeighbors>0,:);
    cliquies = ELSclique(AA);
    node_subGraph= find(nNeighbors>0);
    locs2Group = {};
    nMembers = sum(cliquies);
    while max(nMembers)>1
        [~,idx2Merge] = max(nMembers);
        lMember = logical(cliquies(:,idx2Merge));
        locs2Group{end+1} = node_subGraph(lMember);
        cliquies(lMember,:)=[];
        node_subGraph(lMember)=[];
        
        nMembers = sum(cliquies);
    end
    idxLocsMerged = unique([locs2Group{:}]);
    
    for k = 1:length(fn)
        newLocs.(fn{k}) = locs.(fn{k});
    end
    
    % remove all locs that will be merged
    for k = 1:length(fn)
        newLocs.(fn{k})(idxLocsMerged) = [];
    end
    newLocs.n = ones(length(newLocs.xnm),1);
    for l = 1:length(locs2Group)
        for k = 1:length(fn)
            idxMember = locs2Group{l};
            newLocs.(fn{k})(end+1) = mean(locs.(fn{k})(idxMember));
        end
        newLocs.n(end+1) = length(idxMember);
    end
end