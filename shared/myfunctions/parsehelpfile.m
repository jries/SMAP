function [description,tooltips]=parsehelpfile(filename)
txt=fileread(filename);

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

for k=1:length(indt)-1
    txth=txttt(indt(k)+4:indt(k+1)-1);
    nameend=find(txth<'0' | (txth>'9' & txth<'A'),1,'first'); %allow for a-z, A-Z,0_9,'_'
    hname=txth(1:nameend-1);
    tooltiptxt=reformatline(txth(nameend+1:end));
    if ~strcmp(hname,'Parameters')
        tooltips.(hname)=tooltiptxt;
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