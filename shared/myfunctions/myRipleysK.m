function K2=myRipleysK(x,y,r,box)

% locs(:,1) x-coords
% locs(:,2) y-coords
% r: distance vector: here K(r) values are calculated
% box=[minX maxX minY maxY] here K is calcualted

maxN=5000;
N=length(x);
runs=ceil(N/maxN);
maxNn=floor(N/runs);
maxNn*runs;
% asdf
locs=[x,y];

h=0*r';
% runs=1;
for k=1:runs
    rangek=(k-1)*maxNn+1:k*maxNn;
%     rangek=1:N;
    X1=repmat(locs(rangek,1),1,maxNn);
    Y1=repmat(locs(rangek,2),1,maxNn);
    for l=1:runs
        disp([k,l])

        rangel=(l-1)*maxNn+1:l*maxNn;
        equal=(k==l);
        hx=Rip2(locs(rangek,1:2),locs(rangel,1:2),r,box,equal,X1,Y1);
%         whos hx h
        h=h+hx;
    end
end
h2=h/runs;


K=cumsum(h2);

lambda = maxNn*runs/((box(2)-box(1))*(box(4)-box(3)));
% K = K/lambda/sum(I);
K=K/lambda;
K2=[0; K(1:end-1)];
% K2=K;

% 
% DX = repmat(locs(:,1),1,N)-repmat(locs(:,1)',N,1); %from RipleysK
% DY = repmat(locs(:,2),1,N)-repmat(locs(:,2)',N,1);
% DIST = sqrt(DX.^2+DY.^2);
% % DIST = sort(DIST);
% 
% % DIST(1,:)=[];%corresponds to 2:end. Avoid by subtracting self-distance
% %density
% rbox = min([locs(:,1)'-box(1);box(2)-locs(:,1)';locs(:,2)'-box(3); box(4)-locs(:,2)'] );
% % rbox is distance to box
% I = (rbox>max(r));
% DIST2=DIST(:,I);
% %     if length(I)>0
% %         K(k) = sum(sum(DIST(2:end,I)<dist(k)))/length(I);
% %     end
% 
% 
% h=histc(DIST2(:),r);
% h(1)=h(1)-sum(I); %subtract self-distance =0
% % h(1)=0;
% K=cumsum(h);
% 
% lambda = N/((box(2)-box(1))*(box(4)-box(3)));
% K = K/lambda/sum(I);
% K2=[0; K(1:end-1)];

%same result as RipleysK on random numbers


function h=Rip2(locs1,locs2,r,box,equal,X1,Y1)
% whos locs1 locs2
[N,~] = size(locs2);
% DX = repmat(locs1(:,1),1,N)-repmat(locs2(:,1)',N,1); %from RipleysK
% DY = repmat(locs1(:,2),1,N)-repmat(locs2(:,2)',N,1);
% whos X1
DX = X1-repmat(locs2(:,1)',N,1); %from RipleysK
DY = Y1-repmat(locs2(:,2)',N,1);
DIST = sqrt(DX.^2+DY.^2);
% DIST = sort(DIST);

% DIST(1,:)=[];%corresponds to 2:end. Avoid by subtracting self-distance
%density
rbox1 = min([locs1(:,1)'-box(1);box(2)-locs1(:,1)';locs1(:,2)'-box(3); box(4)-locs1(:,2)'] );
rbox2 = min([locs2(:,1)'-box(1);box(2)-locs2(:,1)';locs2(:,2)'-box(3); box(4)-locs2(:,2)'] );
rbox=min([rbox1; rbox2]);
% whos rbox
% rbox is distance to box
I = (rbox>max(r));
DIST2=DIST(:,I);
%     if length(I)>0
%         K(k) = sum(sum(DIST(2:end,I)<dist(k)))/length(I);
%     end


h=histc(DIST2(:),r);
% h(1)=h(1)-equal*sum(I); %subtract self-distance =0
h=h/(sum(I));
% h(1)=0;
% K=cumsum(h);

% lambda = N/((box(2)-box(1))*(box(4)-box(3)));
% K = K/lambda/sum(I);
% K2=[0; K(1:end-1)];