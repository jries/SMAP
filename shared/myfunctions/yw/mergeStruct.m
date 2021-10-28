function a = mergeStruct(a,b)
    % it requires a and b to have same fields
    fn = fieldnames(a);
    for k = 1:length(fn)
        a.(fn{k}) = [a.(fn{k});b.(fn{k})];
    end
end