function pout=replacefieldsrecursive(pin,strold,strnew)
pout=[];
if isstruct(pin)
    fn=fieldnames(pin);
    for k=1:length(fn)
        pout.(strrep(fn{k},strold,strnew))=replacefieldsrecursive(pin.(fn{k}),strold,strnew);
    end
else
    pout=pin;
end
    