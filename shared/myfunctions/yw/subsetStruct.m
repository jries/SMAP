function a = subsetStruct(a,ind)
    % subsetting the struct having fields with the same length.
    fn = fieldnames(a);
    for k = 1:length(fn)
        a.(fn{k}) = a.(fn{k})(ind);
    end
end