function WriteYamlSimple(sin, file)
    txt=struct2yaml(sin,'')
end

function txt=struct2yaml(p,prefix)
txt="";
if isstruct(p)||isobject(p)
    
    fn=fieldnames(p);
    for k=1:length(fn)
        txt(end+1,1)=string([prefix fn{k} ':']);
        prefixn=[ prefix,'  '];
        to=struct2yaml(p.(fn{k}),prefixn);
        for l=1:length(to)
            txt(end+1,1)=to(l);
        end
    end
else
    if iscell(p)
        to=p{1};
        if isnumeric(to)||islogical(to)
            to=num2str(to);
        end
        for k=2:length(p)
            th=p{k};
            if isnumeric(th)||islogical(th)
                th=num2str(th);
            end
            if ischar(th)
                to=[to ',' th];
            end
        end
        p=to;
    end
    if isnumeric(p)||islogical(p)
        p=num2str(p);
    end
    if ~isempty(p)
        tx=p(1,:);
        if ischar(tx)
            txt=string([prefix ': ' tx]);
        end
    else
        txt=string([prefix ': ']);
    end

end
    txt=txt(2:end);
end