function sout=ReadYamlSimple(file)
ln=readlines(file);
ln=strrep(ln,' - ','  - '); %increase level of arrays
indgood=true(length(ln),1);
for k=1:length(ln)
    ln(k)=removecomment(ln(k));
    if (strip(ln(k),'both'))==""
        indgood(k)=false;
    end
end
sout=recursivesplit(ln(indgood));
end

function strout=recursivesplit(lines)
strout=[];
if isempty(lines)
    return
end
numspace=strlength(lines(1))-strlength(strip(lines(1),'left'));
ind1=1;
ind2=2;
arrayind=1;
while ind1<=length(lines)
    if lines(ind1)==""
        ind1=ind1+1;
        ind2=ind1+1;
        continue
    end
    [name, valout]=parseline(lines(ind1));
    while ind2<=length(lines) && strlength(lines(ind2))-strlength(strip(lines(ind2),'left'))>numspace  
        ind2=ind2+1;
    end
    if ind1+1<=ind2-1  
        if isempty(name)
            linesh=lines(ind1:ind2-1);
            linesh(1)=strrep(linesh(1), '-   - ','  - ');
            val=recursivesplit(linesh);
            strout{arrayind}=val;
            arrayind=arrayind+1;
        else
            val=recursivesplit(lines(ind1+1:ind2-1));
            strout.(name)=val;
        end
        ind1=ind2;
        ind2=ind2+1;
    elseif ~isempty(name)
        strout.(name)=valout;
        ind1=ind1+1;
        ind2=ind1+1;
    else
        %number
        strout{arrayind}=valout;
        arrayind=arrayind+1;
        ind1=ind1+1;
        ind2=ind1+1;      
    end
end
end

function lineout=removecomment(line)
str=char(line);
indcomment=strfind(str,'#'); %remove comments
if isempty(indcomment)
    endline=length(str);
else
    endline=indcomment-1;
end
str=str(1:endline);
lineout=string(str);
end

function [name, valout, novalue,ind,indexlevel]=parseline(line)
indexlevel=1;
novalue=false;
str=char(line);
valout=[];
name=[];
if isempty(str)
    return
end
ind=1;
while(str(ind)==' ')
    ind=ind+1;
end
indcolon=strfind(str,':');
if isempty(indcolon)

    novalue=false;
    inlevel=strfind(str,'- ');
    indexlevel=length(inlevel);
    valout=str2num(str(inlevel(end)+2:end));
    return
end
name=str(ind:indcolon-1);
value=str(indcolon+1:end);
if isempty(value)
    novalue=true;
    valout=[];
    return
end
switch value
    case 'true'
        valout=true;
    case 'false'
        valout=false;
    case 'null'
        valout=[];
    case {'NaN','nan'}
        valout=NaN;
    otherwise
        valout=str2num(value);
        if isempty(valout)
            valout=removejunk(value);
        end
end
end

function strout=removejunk(strin)
strout=strtrim(strin);
if isempty(strout)
    return
end
if strout(1)==''''
    strout(1)=[];
end 
if strout(end)==''''
    strout(end)=[];
end 
end