function out=myimfilter(in,nk)
dn=round((nk-1)/2);
out=myimfilterc(single(in),uint32(dn));

% 
% 
% s=size(in);
% 
% tic
% out=zeros(s(1),s(2));
% for x1=1+dn:s(1)-dn
%     for y1=1+dn:s(2)-dn
%         for x2=x1-dn:x1+dn
%             for y2=y1-dn:y1+dn
%                 out(x1,y1)=out(x1,y1)+in(x2,y2);
%             end
%         end
%     end
% end
% out=out/(2*dn+1)^2;
% t1=toc;
% 
% 
% figure(1)
% subplot(2,2,1)
% imagesc(out)
% colorbar
% title(t1);
% 
% tic

% t2=toc;
% 
% subplot(2,2,2)
% imagesc(out2)
% colorbar
% title(t2);
% subplot(2,2,3)
% imagesc(out2-out)
% colorbar

