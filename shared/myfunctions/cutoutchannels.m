function imout=cutoutchannels(imstack,s)
if nargin<2
    s=readstruct('settings_3D.txt');
end
if isfield(s,'x4pi') %this is a 4pi data set
    numchannel=length(s.x4pi);
    imout=zeros(s.height4pi,s.width4pi*numchannel,size(imstack,3),'like',imstack);
    for k=1:numchannel
%         rangeh=(k-1)*s.height4pi+1:k*s.height4pi;
        rangeh=(k-1)*s.width4pi+1:k*s.width4pi;
        imh=imstack(s.y4pi(k)+1:s.y4pi(k)+s.height4pi,s.x4pi(k)+1:s.x4pi(k)+s.width4pi,:);
        if s.mirror4pi(k)==2 %along x
            imh=imh(end:-1:1,:,:);
        elseif s.mirror4pi(k)==1 %along y
            imh=imh(:,end:-1:1,:);
        end
        imout(:,rangeh,:)=imh;
    end
end
    



