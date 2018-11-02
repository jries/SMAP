function out=prop2struct(in)
if isstruct(in)
    out=in;
    return
end
fn=properties(in);
for k=1:length(fn)
    out.(fn{k})=in.(fn{k});
end
end