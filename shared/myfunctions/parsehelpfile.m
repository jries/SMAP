function [description,tooltips,interpreter]=parsehelpfile(filename)
if exist(filename,'file')
    txt=fileread(filename);
else
    txt=filename;
end
%get interpreter and remove line
indinterpret=strfind(txt,'gui:Interpreter');
if isempty(indinterpret)
    interpreter='tex';
else
   inde=find(txt(indinterpret+16:indinterpret+25)< 'A',1,'first');
   interpreter=txt(indinterpret+16:indinterpret+16+inde-2);
   txt(indinterpret:indinterpret+16+inde-2)=[];
end

indParameters=strfind(txt,'gui:Parameters');
if isempty(indParameters)
    indParameters=strfind(txt,'gui:');
end
description=reformatline(txt(1:indParameters(1)-1));

%tooltips
txttt=txt(indParameters(1):end);
indt=strfind(txttt,'gui:');

%if gui: is used somewhere, remove. Only keep is : is followed by character
indbad=uint8(txttt(indt+4))<65 |uint8(txttt(indt+4))>122;
indt(indbad)=[];
indt(end+1)=length(txttt);
copyf=[];
tooltips=[];
for k=1:length(indt)-1
    txth=txttt(indt(k)+4:indt(k+1)-1);
    nameend=find(txth<'0' | (txth>'9' & txth<'A'),1,'first'); %allow for a-z, A-Z,0_9,'_'
    hname=txth(1:nameend-1);
    tooltiptxt=reformatline(txth(nameend+1:end));
    if txth(nameend)=='=' %copy
        nameendc=find(tooltiptxt<'0' | (tooltiptxt>'9' & tooltiptxt<'A'),1,'first'); 
        if isempty(nameendc)
            copyf.(hname)=tooltiptxt;
        else
            copyf.(hname)=tooltiptxt(1:nameendc-1);
        end
    elseif ~strcmp(hname,'Parameters')
        tooltips.(hname)=tooltiptxt;
    end
end



if ~isempty(copyf) && ~isempty(tooltips)
%     tooltips
    fn=fieldnames(copyf);
    for k=1:length(fn)
        if  isfield(tooltips, copyf.(fn{k}))
            tooltips.(fn{k})=tooltips.(copyf.(fn{k}));
        else 
            display(['check definition for ' fn{k} ' in ' filename])
        end
    end
end
%ADD: gui:h1=h2
end

function out=reformatline(in)
indend=find(in~=newline,1,'last');
out=in(1:indend);
out=strrep(out,[newline newline],'\n');
out(out==newline)=' ';
if ~isempty(out) && out(1)==':'
    out(1)=[];
end
end