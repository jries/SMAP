function savematparts(filen,ld,fields)
saved=trysave(filen,ld);
if saved
    return;
end

% try 
%     ld=ldin.copy;
% catch
%     ld=ldin;
% end
% ld will be destroyed. Copy
if ~iscell(fields)
    fields={fields};
end

% lh=l;
S=struct('type','.','subs',fields);
lh=subsref(ld,S);
% for k=1:length(fields)
%     lh=lh.(fields{k});
% end

[path,file,ext]=fileparts(filen);

fn=fieldnames(lh);
maxsize=2e9; %1gb
ltemp=[];
sizebytes=0;
indfile=1;


partnames={};
for k=1:length(fn)
    vh=lh.(fn{k});
    infh=whos('vh');
    sizebytes=sizebytes+infh.bytes;
    if sizebytes>maxsize
        indfile=indfile+1;
        sizebytes=infh.bytes;
        nameh=[path filesep file '_p' int2str(indfile) ext];
        save(nameh,'ltemp');
        partnames{indfile-1}=nameh;
        fnsaved=fieldnames(ltemp);
        for l=1:length(fnsaved)
            lh.(fnsaved{l})=[];
        end
        
        ltemp=[];
    end
    ltemp.(fn{k})= vh;
    

end
lds=subsasgn(ld,S,lh);

saved=trysave([path filesep file ext],struct('lds',lds,'S',S,'partnames',{partnames}));
if ~saved
    indfile=indfile+1;
    nameh=[path filesep file '_p' int2str(indfile) ext];
    save(nameh,'ltemp');
    partnames{indfile-1}=nameh;
    lds=subsasgn(ld,S,[]); 
    saved=trysave([path filesep file ext],struct('lds',lds,'S',S,'partnames',{partnames}));
    if ~saved
        save([path filesep file ext],'lds','S','partnames','-v7.3')
    end
end
%try to save all, if it doesnt work, save ltemp and rest

end


function saved=trysave(file,ls)
saved=false;
try
    save(file,'-struct','ls','-v7');    
    [msg,msgid]=lastwarn;
    if strcmp(msgid,'MATLAB:save:sizeTooBigForMATFile')
        saved=false;
        lastwarn('cleared');
    else
        saved=true;
    end
end
end