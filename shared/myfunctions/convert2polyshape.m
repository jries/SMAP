function p=convert2polyshape(roiin, lw)
    if isempty(roiin)
        p=[];
    elseif isstruct(roiin)
        p=rect2poly(roiin.getPosition*1000);
    elseif isnumeric(roiin)&& numel(roiin) ==4 %rectangle
        p=rect2poly(roiin);
    elseif isa(roiin,'imline')
        coord=roiin.getPosition*1000;
        
        vp=coord(2,:)-coord(1,:);
         angle=atan2(vp(2),vp(1));
        len=sqrt(sum(vp.^2));

        p=rect2poly([0 0 len/2 lw/2]);
        [xs,ys]=rotcoord(p.Vertices(:,1),p.Vertices(:,2), -angle);
        mp=mean(coord,1); 
        p=polyshape(mp(1)+xs,mp(2)+ys);

        

%         vs=[vp(2), -vp(1)];
%         vs=vs./abs(vs);
%         x=[coord(1,1)+ lw* vs(1);coord(2,1)+ lw* vs(1);coord(2,1)- lw* vs(1);coord(1,1)- lw* vs(1)];
%         y=[coord(1,2)+ lw* vs(2);coord(2,2)+ lw* vs(2);coord(2,2)- lw* vs(2);coord(1,2)- lw* vs(2)];
%     
    elseif isa(roiin,'imellipse')
        coord=roiin.getPosition;
         Polygons=Ellipses2Polygons([coord([3 4 1 2]) 0]*1000,1000);
         p=polyshape(Polygons{1});

     elseif isa(roiin,'imrect')
         p=rect2poly(roiin.getPosition*1000);
    elseif isa(roiin,'polyshape')
        p=roiin;
%         coord=roih.getPosition*1000;
%         p=polyshape(coord(:,1),coord(:,2));
    else
        roiin
        disp('roi not recognized')
    end

end

function p=rect2poly(c)
coord=[c(1)-c(3), c(2)-c(4);c(1)+c(3), c(2)-c(4);c(1)+c(3), c(2)+c(4);c(1)-c(3), c(2)+c(4)];
p=polyshape(coord(:,1),coord(:,2));
end