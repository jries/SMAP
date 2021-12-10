function txt=WriteYamlSimple(file,sin)
    txt=struct2yaml(sin,'');
    txto=txt(2:end);
fid = fopen(file, 'w');
fprintf(fid, '%s\n', txto);
fclose(fid);
end

function txt=struct2yaml(p,prefix)
txt="";
if isstruct(p)||isobject(p)
    
    fn=fieldnames(p);
    for k=1:length(fn)
        txt(end+1,1)=string([prefix fn{k} ':']);
        prefixn=[ prefix,'  '];
        to=struct2yaml(p.(fn{k}),prefixn);
        txt(end,1)=txt(end,1)+" "+to(1);
        for l=2:length(to)
            txt(end+1,1)=to(l);
        end
    end
else
    if iscell(p) || (isnumeric(p) && length(p)>1)
        s=cells2string(p, [prefix(3:end) '  ']);
        txt(end+1:end+length(s))=s;
%         to=p{1};
%         if isnumeric(to)||islogical(to)
%             to=num2str(to);
%         end
%         for k=2:length(p)
%             th=p{k};
%             if isnumeric(th)||islogical(th)
%                 th=num2str(th);
%             end
%             if ischar(th)
%                 to=[to ',' th];
%             end
%         end
%         p=to;
    else
        txt(1,1)=makestring(p);
    end

end

end

function s=cells2string(p, prefix)
s="";
for k=1:length(p)
    if iscell(p)
        pk=p{k};
    else
        pk=p(k);
    end
%     prefixh=prefix;
%     if k>1
%         prefixh(3)=' ';
%     end
    prefix1=prefix;prefix1(end-1)='-';
    if iscell(pk)
%         sh=cells2string(p{k},[prefixh '- ']);
        sh=cells2string(pk,'  ');
        
        s(end+1,1)=string(prefix1)+sh(1);
        for ks=2:length(sh)
            s(end+1,1)=string(prefix)+sh(ks);
        end
    else
        ps=makestring(pk);
        if isempty(ps)
            ps="null";
        end
        s(end+1,1)= string(prefix1)+ ps;
    end

end
s=s(2:end,1);
end

function sout=makestring(in)
if isempty(in)
    sout="null";
elseif isfloat(in)
    if round(in) ~= in
        sout=string(num2str(in,'%2.1f'));
    else
        sout=string(in);
    end
elseif ischar(in)
    sout=string(in);
else
%     in
    sout=string(in);
end

end