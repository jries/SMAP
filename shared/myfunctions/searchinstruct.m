function outstr=searchinstruct(instr,fieldname,compfun)
if nargin<3
    compfun=@strcmp;
end
outstr='';
if ~isstruct(instr)
    return
end
fn=fieldnames(instr);
ind=find(compfun(fn,fieldname));
if ~isempty(ind)
    outstr=(fn{ind});
    return
end

for k=1:length(fn)
    sstr=searchinstruct(instr.(fn{k}),fieldname,compfun);
    if ~isempty(sstr)
        outstr=[fn{k} '.' sstr];
    end
end
