function img=insertLineImage(xin,yin,img,val)
if nargin<4
    val=1;
end
for k=2:length(xin)
    x2=xin(k);x1=xin(k-1);
    y2=yin(k);y1=yin(k-1);

    % Bresenham's line algorithm
    dx = abs(x2 - x1);
    dy = abs(y2 - y1);
    sx = sign(x2 - x1);
    sy = sign(y2 - y1);
    
    err = dx - dy;
    
    x = x1;
    y = y1;
    
    while true
        % Set the pixel value to 1 (or any desired value)
        img(round(y), round(x)) = val;
        
        % Check if the end point is reached
        if x == x2 && y == y2
            break;
        end
        
        % Update error and coordinates
        e2 = 2 * err;
        if e2 > -dy
            err = err - dy;
            x = x + sx;
        end
        if e2 < dx
            err = err + dx;
            y = y + sy;
        end
    end
end