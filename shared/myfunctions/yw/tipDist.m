function dist = tipDist(tips)
    x1 = tips(1,1);
    x2 = tips(2,1);
    y1 = tips(1,2);
    y2 = tips(2,2);
    
    dist = sqrt((y2-y1)^2+(x2-x1)^2);
end