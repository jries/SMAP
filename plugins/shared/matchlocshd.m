function [iA,iB,uiA,uiB]=matchlocshd(xA,yA,xB,yB,d,maxd)
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
% matched=false(size(ixBs));
% xAs==xA(ixAs)

% iA=zeros(lA,1);iB=zeros(lB,1);
dmatrix=zeros(lA,lB);

for k=1:lA
        l=lold;
        while  l<=lB && xBs(l)+d(1)<xAs(k)-maxd% move till possible hit
            l=l+1;
            if l>lB %no match
                break
            end
        end
        lold=l; %
%         dold2=maxd^2;
%         found=0;
        while l<=lB && xBs(l)+d(1)<xAs(k)+maxd  %all possible hits
            dh=((xBs(l)+d(1)-xAs(k))^2+(yA(ixAs(k))-yB(ixBs(l))-d(2))^2);
            if dh<maxd^2
                dmatrix(k,l)=1/dh;
            end
%             if dnew2<dold2 &&~matched(l)
%                 iA(rx)=ixAs(k);
%                 iB(rx)=ixBs(l);
%                 matched(l)=true; %associate only onces
%                 dold2=dnew2;
%                 found=1;
%             end
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
%         if found==1
%             rx=rx+1;
%         end
end

dmatrixs=sparse(dmatrix);
[dm,sortindd]=sort(dmatrixs(:),'descend');

% numcandidates=find(full(dm)==0,1,'first');
numcandidates=nnz(dmatrixs);
[ik,il]=ind2sub(size(dmatrix),sortindd(1:numcandidates));
ig=1;

iAi=zeros(numcandidates,1);iBi=zeros(numcandidates,1);
% for k=1:numcandidates
%     if dmatrix(ik(k),il(k))>0
%         iAi(ig)=ik(k);
%         iBi(ig)=il(k);
%         dmatrix(ik(k),:)=0;
%         dmatrix(:,il(k))=0;
%         ig=ig+1;
%     end
%     
% end
ikm=true(size(dmatrix,1));
ilm=true(size(dmatrix,2));
for k=1:numcandidates
    if ikm(ik(k)) && ilm(il(k))
        iAi(ig)=ik(k);
        iBi(ig)=il(k);
        ikm(ik(k)) =0;
        ilm(il(k))=0;
%         ik(ik==ik(k))=0;
%         il(il==il(k))=0;
        ig=ig+1;
    end
    
end


% [vm,indm]=max(dmatrix,[],2);
% numfound=sum(vm>0);
% iAi=zeros(numfound,1);iBi=zeros(numfound,1);
% ig=1;
% [~,vsind]=sort(vm,'descend');
% for k=1:lA
%     indA=vsind(k);
%     indB=indm(vsind(k));
%     if vm(indA)>0 && ~any(iBi==indB)
%         iAi(ig)=(indA);iBi(ig)=(indB);
%         ig=ig+1;
%     end
% end
% 


iAi(ig:end)=[];
iBi(ig:end)=[];
% 
% figure(88);
% plot(vertcat(xAs(iAi)',xBs(iBi)'),vertcat(yA(ixAs(iAi))',yB(ixBs(iBi))'))

iA=ixAs(iAi);
iB=ixBs(iBi);
uiA=1:lA;
uiB=1:lB;
uiA(iA)=[];
uiB(iB)=[];
if 0
figure(89);
hold off
plot(vertcat(xA(iA)',xB(iB)'),vertcat(yA(iA)',yB(iB)'),'r-', xA(uiA),yA(uiA),'go', xB(uiB),yB(uiB),'gd')
hold on
plot(vertcat(xA(1:end-1)',xB'),vertcat(yA(1:end-1)',yB'),'b-')
end
end

            