function destination=copystructReduce(source,ind,fn)
if ~isempty(source)
    if nargin<3
        fn=fieldnames(source);
    end
    for k=1:length(fn)
        if length(ind)==length(source.(fn{k}))||~islogical(ind)
            destination.(fn{k})=source.(fn{k})(ind);
        else
            destination.(fn{k})=source.(fn{k});
        end
    end
end
