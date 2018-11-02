function out=fastisprop(object,field)
    fn=fieldnames(object);
    out= sum(strcmp(fn,field));
end