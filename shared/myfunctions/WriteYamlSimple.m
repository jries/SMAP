function txt=WriteYamlSimple(sin, file)
    txt=struct2yaml(sin,'');
    txto=txt(2:end);
fid = fopen(file, 'wt');
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
    if iscell(p)
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
        txt(1,1)=string(p);
    end
%     if isnumeric(p)||islogical(p)
%         p=num2str(p);
%     end
%     if ~isempty(p)
%         tx=p(1,:);
%         if ischar(tx)
%             txt=string([prefix ': ' tx]);
%         end
%     else
%         txt=string([prefix ': ']);
%     end

end
%     txt=txt(2:end);
end

function s=cells2string(p, prefix)
s="";
for k=1:length(p)
%     prefixh=prefix;
%     if k>1
%         prefixh(3)=' ';
%     end
prefix1=prefix;prefix1(end-1)='-';
    if iscell(p{k})
%         sh=cells2string(p{k},[prefixh '- ']);
        sh=cells2string(p{k},'  ');
        
        s(end+1,1)=string(prefix1)+sh(1);
        for ks=2:length(sh)
            s(end+1,1)=string(prefix)+sh(ks);
        end
    else
        ph=string(p{k});
        if isempty(ph)
            ph="null";
        end
        s(end+1,1)= string(prefix1)+ ph;
    end

end
s=s(2:end,1);
end