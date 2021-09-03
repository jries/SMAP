function img3d = NPCImageModel(locs,par)
    if isempty(locs)
        img3d = [];
        img3d.image = [];
    end
    p.size = 300;
    p.r = 105/2;
    p.plotSlice = 1;
    p.std = [10 10 12];
    
    [x,y] = meshgrid(-(p.size/2):(p.size/2), -(p.size/2):(p.size/2));
    v = simpleRing(x(:),y(:),p.r, 10);

    [x,y] = meshgrid(1:p.size+1, 1:p.size+1);
    twoDInd = sub2ind([p.size+1 p.size+1],y(:),x(:));
    twoD = zeros([p.size+1 p.size+1]);

    twoD(twoDInd) = v;
    % twoD = filter2(hl, twoD);
    figure;imagesc(twoD)
    threeD = zeros([p.size+1 p.size+1 p.size+1]);
    threeD(:, :, p.size/2-5:p.size/2+5) = repmat(twoD, [1 1 11]);
    
%     img3d = imgaussfilt3(threeD,p.std);
    img3d = gauss3filter(threeD,p.std);
    img3d = double(img3d);
    %figure; isosurface(threeD)
    if p.plotSlice
        figure; slice(img3d, 150,150,150)
    end
    
end

function v = simpleRing(x,y,r,ran)
        d = x.^2 + y.^2;
        v = d < r.^2+ran^2 & d > r.^2-ran^2;
end