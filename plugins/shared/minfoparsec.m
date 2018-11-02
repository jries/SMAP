function info=minfoparsec(minfo,ctab)
if isempty(minfo)
    info=[];
    return
end
% try
c=ctab;
info=setdefault(c);
%determine which chip
c.Properties.VariableNames{1}='tag';
ind=find(c.tag==1);
indnames=[];
for k=1:length(ind)
    [cam,indcamid]=getvald(minfo,c.camId{ind(k)});
    if ~isempty(cam)
    info.camId=cam;
    indnames=ind(k);
    break
    end
end

if ~isempty(indnames) %cam found

%determine MMCamname
k=indcamid;
while minfo(k)~='"'
    k=k-1;
end
MMCamname=minfo(k+1:indcamid-1);

%find parameters in minfo
info.port=getvald(minfo,c.port{indnames},MMCamname);
info.preAmp=getvald(minfo,c.preAmp{indnames},MMCamname);
info.readoutrate=getvald(minfo,c.readoutrate{indnames},MMCamname);
info.exposure=str2double(getvald(minfo,c.exposure{indnames},MMCamname));
info.emgain=str2double(getvald(minfo,c.emgain{indnames},MMCamname));
if isnan(info.emgain)
    info.emgain=1;
end
roi=str2num(getvaldroi(minfo,c.roi{indnames}));
info.roi=roi(:)';

comment=getvald(minfo,'Comment');
if ~isempty(comment)
    info.comment=comment;
end
% find correct line
indcal=-1;
camIds=c{:,2};
for k=1:length(camIds)
    if streq(cam,camIds{k})&&streq(info.port,c.port{k,1})&&streq(info.preAmp,c.preAmp{k,1})...
            &&streq(info.readoutrate,c.readoutrate{k,1})
        indcal=k;
    end
end
if indcal==-1
    disp('these settings have not been calibrated')
else
%associate last values
info.conversion=str2double(c.conversion{indcal,1});
info.offset=str2double(c.offset{indcal,1});
info.temperature=str2double(c.temperature{indcal,1});
info.cam_pixelsize_um=str2double(c.pixsize{indcal,1});
end
% 
% else  %default parameters
%     info=setdefault(c);
end


%get real time interval
[f1,t1]=minfotime(minfo,1);
[f2,t2]=minfotime(minfo,-1);
if ~isempty(t1)&&~isempty(t2)
timeinterval=(t2-t1)/(f2-f1);
else
    timeinterval=1;
end
if ~isinf(timeinterval)
info.timediff=timeinterval;
end


%used names, remove later
info.gain=info.preAmp;
if isfield(info,'camId')
info.chip=info.camId;
end

comment=getvald(minfo,'Comment');
if ~isempty(comment)
    info.comments=comment;
end
% info.comments=getvald(minfo,'Comment','');

%PI position
if isfield(c,'PiezoZ')
    info.PIz=str2double(getvald(minfo,c.PiezoZ{indnames},''));
% else
%     info.PIz=0;
end
% catch
%     disp('metadata could not be parsed')
% end

function info=setdefault(c)
%     disp('unknown camera')
    indnames=find(strcmp(c.camId,'Default'));
    info.port=c.port{indnames};
    info.preAmp=c.preAmp{indnames};
    info.readoutrate=c.readoutrate{indnames};
    info.exposure=str2double(c.exposure{indnames});
    info.emgain=str2double(c.emgain{indnames});
    info.roi=c.roi{indnames};
%     info.roi=[0 0 0 0];
    info.conversion=str2double(c.conversion{indnames,1});
    info.offset=str2double(c.offset{indnames,1});
    info.temperature=str2double(c.temperature{indnames,1});
    info.cam_pixelsize_um=str2double(c.pixsize{indnames,1});


function [f,t]=minfotime(minfo,direction)
if direction==-1
    indf=strfindfast(minfo,'ElapsedTime-ms',1,direction);
    indf=strfindfast(minfo,'FrameKey',length(minfo)-indf,direction);
else
indf=strfindfast(minfo,'FrameKey',1,direction);
end
if isempty(indf)
    t=[];
    f=[];
    return
end
if direction==-1
    indf=indf-8;
end
sf=minfo(indf:indf+15);

indminus=strfind(sf,'-');
numstr=sf(indminus(1)+1:indminus(2)-1);
f=str2double(numstr);
indt=strfindfast(minfo,'ElapsedTime-ms',indf);
t1s=getval(minfo,indt);
t=str2double(t1s(3:end-1));

function txt=getval(minfo,index,sep1,sep2)
if isempty(index)
    txt='';
    return
end
if nargin==2
    sep1='"';
    sep2='"';
end

while index<length(minfo)&& minfo(index)~=sep1
    index=index+1;
end
i1=index+1;
index=index+1;
while index<length(minfo) && minfo(index)~=sep2
    index=index+1;
end
txt=minfo(i1:index-1);

function [txt,ind]=getvald(minfo,finds,camname,startind)
if nargin==2
    camname=[];
end
if nargin<4
    startind=1;
end
searchstr=[camname  finds '":'];
ind=strfindfast(minfo,searchstr,startind);
txt=getval(minfo,ind+length(searchstr));

function [txt,ind]=getvaldroi(minfo,finds,camname,startind)
if nargin==2
    camname=[];
end
if nargin<4
    startind=1;
end
searchstr=[camname  finds '":'];
% ind=strfindfast(minfo,searchstr,startind);
% txt=getval(minfo,ind+length(searchstr));%,'[',']');
% txt=strrep(txt,'-',' ');
    ind=strfindfast(minfo,searchstr,startind);
    txt1=getval(minfo,ind+length(searchstr));%,'[',']');
    txt1=strrep(txt1,'-',' ');
   
   if isempty(str2num(txt1))
       txt2=getval(minfo,ind+length(searchstr),'[',']');
       txt=txt2;
   else
       txt=txt1;
   end



function out=streq(a,b)
out=~isempty(strfind(a,b));
