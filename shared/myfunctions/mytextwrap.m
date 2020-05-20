function textout= mytextwrap(textin, width, margin,newlinechar)
if nargin<3
    margin=5;
end
if nargin<4
    newlinechar='\n';
end
textout='';
ind1=1;
ind2=0;
while ind1<length(textin)
ind2=min(length(textin),ind1+width);
texth=textin(ind1:ind2);
indf1=strfind(texth,newline);
indf2=strfind(texth,'\n');
if ~isempty(indf1)
    textout=[textout textin(ind1:ind1+indf1-2) newlinechar];
    ind1=ind1+indf1;
elseif ~isempty(indf2)
    textout=[textout textin(ind1:ind1+indf2-2) newlinechar];
    ind1=ind1+indf2+1;
else
    indsp=strfind(texth,' ');
    if ~isempty(indsp) && width-indsp(end)<=margin
        indend=ind1+indsp(end)-2;
        textout=[textout textin(ind1:indend) newlinechar];
        ind1=indend+1;
    else
        textout=[textout textin(ind1:ind2) newlinechar];
        ind1=ind2+1;
    end
end
end
if length(textout)>2
textout(end-1:end)=[];
end

