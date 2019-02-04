function [iA,iB,uiA,uiB]=matchlocs(xA,yA,xB,yB,d,maxd)
%only one same frame
%x=xA(iA)=xB(iB)
%d=[dx dy]
%maxd: maximum allowed difference for match
% uiA,uiB:unmatched
if isempty(d)
    d=[0,0];
end
lA=length(xA); %xA,yA: same length
lB=length(xB);


rx=1;
lold=1;
[xAs,ixAs]=sort(xA);
[xBs,ixBs]=sort(xB);
matched=false(size(ixBs));
% xAs==xA(ixAs)

iA=zeros(lA,1);iB=zeros(lB,1);

for k=1:lA
        l=lold;
        while  l<=lB&&xBs(l)+d(1)<xAs(k)-maxd% move till possible hit
            l=l+1;
            if l>lB %no match
                break
            end
        end;
        lold=l; %
        dold2=maxd^2;
        found=0;
        while l<=lB&&xBs(l)+d(1)<xAs(k)+maxd  %all possible hits
            dnew2=(xBs(l)+d(1)-xAs(k))^2+(yA(ixAs(k))-yB(ixBs(l))-d(2))^2;
            if dnew2<dold2 &&~matched(l)
                iA(rx)=ixAs(k);
                iB(rx)=ixBs(l);
                matched(l)=true; %associate only onces
                dold2=dnew2;
                found=1;
            end
            l=l+1;
            %idea:
            %per frame: matrix all agianst all
            % if two are closer than distance: put distance in matrix
            % associate lowest entry
            % remove, associate next
            % or to speed up: when matched: store d2 and both indices
            % 
           
%             if abs(yA(ixAs(k))-yB(ixBs(l))+d(2))<maxd %true hit
%                 %do something here
%                 iA(rx)=ixAs(k);
%                 iB(rx)=ixBs(l);
%                 rx=rx+1;
%                 break %continure with next particle A
%             else
%                 l=l+1;
%             end
        end
        if found==1
            rx=rx+1;
        end
end
iA(rx:end)=[];
iB(rx:end)=[];


uiA=1:lA;
uiB=1:lB;
uiA(iA)=[];
uiB(iB)=[];
            