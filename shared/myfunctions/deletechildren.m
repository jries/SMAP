function deletechildren(obj)
if isprop(obj,'children')
    if ~isempty(obj.children)
    fn=fieldnames(obj.children);
    for k=1:length(fn)
        deletechildren(obj.children.(fn{k}));
    end
    end
end
if isprop(obj,'handle')
    delete(obj.handle)
end
    delete(obj);

end