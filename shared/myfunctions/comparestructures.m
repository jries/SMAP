function comparestructures(a,b)
fna=fieldnames(a);
fnb=fieldnames(b);
string(setdiff(fna,fnb)) 
string(setdiff(fnb,fna))
fn=intersect(fna,fnb);
for k=1:length(fn)
    if isstruct(a.(fn{k}))
        comparestructures(a.(fn{k}),b.(fn{k}));
    elseif isnumeric(a.(fn{k}))
        if a.(fn{k}) ~= b.(fn{k})
            disp(fn{k},a.(fn{k}),b.(fn{k}))
        end
    else 
        a.(fn{k}),b.(fn{k})
    end
end