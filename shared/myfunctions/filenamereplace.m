function fout=filenamereplace(fin,fref,ftar)

% provide matching pair fref, ftar
% provide fin, then fout will be matched to ftar

% directories: if different: replace old by new
%files: look at parts between _ and .

fref=['/' strrep(fref,filesep,'/')];ftar=['/' strrep(ftar,filesep,'/')];fin=[strrep(fin,filesep,'/')];
frefind=strfind(fref,'/');frefind=[frefind];
ftarind=strfind(ftar,'/');ftarind=[ftarind];
finind=strfind(fin,'/');finind=[finind];
% look at file
fileref=(fref(frefind(end)+1:end));
filetar=(ftar(ftarind(end)+1:end));
filein=(fin(finind(end)+1:end));
fileout=replacename(filein,fileref,filetar);


% look at path
pathout=fin(1:finind(end));
    for k=0:min([length(frefind) length(ftarind) length(finind)])-2
        pr=fref(frefind(end-k-1)+1:frefind(end-k)-1);
        pt=ftar(ftarind(end-k-1)+1:ftarind(end-k)-1);
        if ~strcmp(pr,pt)
            ph=pathout(finind(end-k-1)+1:finind(end-k)-1);
            pht=replacename(ph,pr,pt);
            pathout=[pathout(1:finind(end-k-1)) pht pathout(finind(end-k):end)];
%             pathout(finind(end-k-1)+1:finind(end-k)-1)=pt;
        end
    end
fout=[pathout fileout];
% fout=pathout;
end

function fileout=replacename(filein,fileref,filetar)
    fileout='';
if ~strcmp(fileref,filetar)
    % search for _ and .
    frr=fileref;frt=filetar;fri=filein;
%     frr=['_' strrep(fileref,'.','_') '_'];frt=['_' strrep(filetar,'.','_') '_'];fri=['_' strrep(filein,'.','_') '_'];
    frri=[1 strfind(frr,"_"|".") length(frr)+1];frti=[1 strfind(frt,"_"|".") length(frt)+1];frii=[1 strfind(fri,"_"|".") length(fri)+1];
    for k=0:min([length(frri) length(frti) length(frii)])-2
        pr=frr(frri(end-k-1):frri(end-k)-1);pt=frt(frti(end-k-1):frti(end-k)-1);pi=fri(frii(end-k-1):frii(end-k)-1);
%         pi=fileout(frii(end-k-1):frii(end-k)-2);
% if the same: use this part
%if different: use the corresponding part from filein
        if strcmp(pi,pr)
            if strcmp(pi,pt)
                %keep either
                fileout=[pi fileout];
            else
                fileout=[pt fileout];
%                 disp('pattern do not seem to match. Output reference file.')
%                 fileout=fileref;
%                 return
                %warning
            end
            
        else
            if strcmp(pr,pt)
                fileout=[pi fileout];
                %keep in
            else
                %keep tar
                fileout=[pt fileout];
            end            
%             fileout=[fileout(1:frii(end-k-1)-1) pt fileout(frii(end-k)-1:end)];
%             fileout(frii(end-k-1):frii(end-k)-2)=pt;
        end
    end
    %append missing target parts
    
end
end