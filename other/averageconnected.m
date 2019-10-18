n=zeros(64);
sim=size(n);
cutoff=0.1;
npos=1500;
linind=ceil(rand(npos,1)*size(n,1)*size(n,2));
n(linind)=rand(npos,1);

n(n<cutoff)=0;
% n=n>cutoff;

nstart=n;
normcounter=ones(size(n));

dxmap=rand(size(n));
dymap=rand(size(n));
dzmap=rand(size(n));
dmaxlat=.8;
dmaxax=.8;
dxmapc=dxmap;

weights=rand(size(n)); %for weighted average. P or photons or CRLB (then for each map individually)
%can we just multiply all maps with weights? for normalization we have to
%add up weights.
dxmapc=dxmapc.*weights;
tic
for x=2:sim(1)-1
    for y=2:sim(2)-1
        if n(x,y)<cutoff
            continue
        end
        for nx=-1:1
            for ny=-1:1
                if n(x+nx, y+ny)<cutoff || (nx==0 &&ny==0)
                    continue
                end
                if n(x+nx, y+ny) > n(x,y)
                    if (dxmap(x+nx, y+ny)-dxmap(x, y))^2+(dymap(x+nx, y+ny)-dymap(x, y))^2<dmaxlat^2 && (dzmap(x+nx, y+ny)-dzmap(x, y))^2<dmaxax^2
                        
                        dxmapc(x+nx,y+ny)=dxmapc(x+nx,y+ny)+dxmapc(x,y); %add map x,y to map x+xn,y+yn %here: add weighting!
                        
                        normcounter(x+nx,y+ny)=normcounter(x+nx,y+ny)+1;
                        %with weights:
                        weights(x+nx,y+ny)=weights(x+nx,y+ny)+weights(x,y);
                        n(x,y)=0; %remove
                    end
                else
                    if (dxmap(x+nx, y+ny)-dxmap(x, y))^2+(dymap(x+nx, y+ny)-dymap(x, y))^2<dmaxlat^2 && (dzmap(x+nx, y+ny)-dzmap(x, y))^2<dmaxax^2
                        dxmapc(x,y)=dxmapc(x,y)+dxmapc(x+nx,y+ny);
                        weights(x,y)=weights(x+nx,y+ny)+weights(x,y);
                        normcounter(x,y)=normcounter(x,y)+1;
                        n(x+nx,y+ny)=0;
                    end
                end
            end
        end
    end
end
toc
% dxmapfinal=dxmapc./normcounter;
dxmapfinal=dxmapc./weights;
  figure(88);
imagesc(vertcat(horzcat(nstart,n,weights),horzcat(dxmap,dxmapc,dxmapfinal)))

%use normcounter or n to find positions and return list
