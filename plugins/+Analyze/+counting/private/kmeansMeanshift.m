function posnew=kmeansMeanshift(x,y,sigma)
X=[x,y];
C=kmeanssplit(x,y);
Cms=meanshift(X,C,sigma,sigma/50);
posnew=spacedApart(Cms,sigma/3);

if 0
figure(25)
plot(x,y,'b.')
hold on
plot(C(:,1),C(:,2),'ro')
plot(Cms(:,1),Cms(:,2),'r+')
plot(posnew(:,1),posnew(:,2),'kx')
hold off
end


function Cms=meanshift(X,C,sigma,cutoff)
sC=size(C);


Cms=zeros(sC);
for k=1:sC(1)
    meanold=C(k,:);
%     inrange=((X(:,1)-C(k,1)).^2+(X(:,2)-C(k,2)).^2)<sigma^2*4;
%     X2=X(inrange,:);
% X2=X;
    while 1
    meannew=wmean(X,meanold,sigma);
        if sum((meannew-meanold).^2)<cutoff^2
            break
        else
            meanold=meannew;
        end
    end
    Cms(k,:)=meannew;
end

function posnew=spacedApart(posold,cutoff)
% ind=true(length(posold),1);
ind=zeros(length(posold),1);
clustern=1;
for k=1:length(posold)
%     d2=((posold(k,1)-posnew(:,1)).^2)+((posold(k,2)-posnew(:,2)).^2);
    
    if ind(k)==0
        d2=((posold(k,1)-posold(:,1)).^2)+((posold(k,2)-posold(:,2)).^2);
        
        ind(d2<cutoff^2)=clustern;
        clustern=clustern+1;
    end
    
end
posnew=zeros(clustern-1,2);
for k=1:clustern-1
    indk=ind==k;
    if sum(indk)==1
        posnew(k,:)=posold(indk,:);
    else
    posnew(k,:)=mean(posold(indk,:),1);
    end
end

%group results close together



function meannew=wmean(X,meanold,sigma)

Xk(:,2)=X(:,2)-meanold(2);
Xk(:,1)=X(:,1)-meanold(1);

K=kernel(Xk,sigma);
meannew(1)=sum(X(:,1).*K)/sum(K);
meannew(2)=sum(X(:,2).*K)/sum(K);

function out=kernel(X,sigma)
s2=sigma^2;
out=exp((-X(:,1).^2-X(:,2).^2)/2/s2);

function C=kmeanssplit(x,y)
chunk=5e2;
numberlocs=length(x);
numbertiles=ceil(sqrt(numberlocs/chunk))
dx=(max(x)-min(x))/numbertiles;dy=(max(y)-min(y))/numbertiles;
mx=min(x);my=min(y);

C=[];
tic
for k=1:numbertiles
    disp(k)
    for l=1:numbertiles
        igb=x>mx+(k-1)*dx&x<mx+k*dx&y>my+(l-1)*dy&y<mx+l*dy;
        np=round(sum(igb)/6);
        
        if np>0
        X=[x(igb),y(igb)];
        
        [~,C1]=kmeans(X,np,'MaxIter',25);
        
        C=[C;C1];
        end
    end
end
toc
