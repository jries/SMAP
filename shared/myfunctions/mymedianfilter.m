function out= mymedianfilter(in,rho)
% tic
%median filter using a rho x rho square for dim 1,2 and all of dim 3
dr=floor((rho-1)/2);
s=size(in);
if length(s)==2;
    s(3)=1;
end

    

lz=(2*dr+1)^2*s(3);
xs=dr+1:2*dr+1:s(1)-dr;
ys=dr+1:2*dr+1:s(2)-dr;

outim=zeros(length(xs),length(ys),class(in));
if s(3)>0
reorderim=zeros(lz,length(ys),class(in));
for x=1:length(xs)
%     indy=1;
    for y= 1:length(ys) 
        si=in(xs(x)-dr:xs(x)+dr,ys(y)-dr:ys(y)+dr,:);
        reorderim(:,y)=si(:);
%         indy=indy+1;
    end
    outim(x,:)=fast_median(reorderim);

end
end
[X,Y]=meshgrid(xs,ys);
[XI,YI]=meshgrid(1:s(1),1:s(2));

out=interp2(X,Y,outim',XI,YI)';


% out(isnan(out))=max(out(:));
% for k=1:dr
% out(k,:)=out(dr+1,:);
% out(end-k+1,:)=out(end-dr,:);
% out(:,k)=out(:,dr+1);
% out(:,end-k+1)=out(:,end-dr);
% end
% size(out)
% toc
% figure(1);imagesc(out)
%       figure(2);imagesc(outim)