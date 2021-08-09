function subFiles = fileSubsets(files, exp, which2Check, exludsive)
    % find out the files with the specified rep number  
    lKept = regexp({files.(which2Check)}, exp);
    if exludsive
    	lKept = cellfun(@isempty,lKept);
    else
        lKept = ~cellfun(@isempty,lKept);
    end
    subFiles  = files(lKept);
end