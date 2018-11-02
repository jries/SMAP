function [beadnum,numlocs,dist,numlocsred,inreduced]=associatenearest(mx,my,pos,maxd)
%associate only to nearest position
if nargin==3
    maxd2=inf;
else
    maxd2=maxd^2;
end

if ~isfield(pos,'x')
    posh.x=pos(:,2);
    posh.y=pos(:,3);
else
    posh=pos;
end

inreduced=false(length(posh.x),1);
beadnum=zeros(length(posh.x),1);
dist=beadnum;
numlocs=zeros(length(mx),1);
numlocsred=zeros(length(mx),1);
for k=1:length(posh.x)
    
    [d2,minind]=min((posh.x(k)-mx).^2+(posh.y(k)-my).^2);
    beadnum(k)=minind;
    numlocs(minind)=numlocs(minind)+1;
    if d2<maxd2
    numlocsred(minind)=numlocsred(minind)+1;
    inreduced(k)=true;
    end
    dist(k)=sqrt(d2);
end

% 
% beadnum=zeros(length(posh.x),1);
% numlocs=beadnum;
% distance2=beadnum+100000000;
% for k=1:length(mx)
%     r2=(posh.x-mx(k)).^2+(posh.y-my(k)).^2;
%     ind=r2<maxd^2&r2<distance2;
%     distance2(ind)=r2(ind);
%     
%     beadnum(ind)=k;
%     numlocs(k)=sum(ind);
% end