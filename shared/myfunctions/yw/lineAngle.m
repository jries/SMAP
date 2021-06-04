function degree = lineAngle(tips)
    x1 = tips(1,1);
    x2 = tips(2,1);
    y1 = tips(1,2);
    y2 = tips(2,2);
    
    rad = atan((y2-y1)/(x2-x1));
    degree = rad2deg(rad);
end