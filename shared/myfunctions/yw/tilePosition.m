function [r,c] = tilePosition(k, dim, dir)
% convert the index k into subsrciption given the direction of increase
% 'dir' and the dimension 'dim'
    switch dir
        case 'down'
            sitePerUnit = dim(1);
        case 'right'
            sitePerUnit = dim(2);
    end
    i = ceil(k/sitePerUnit);
    j = rem(k,sitePerUnit);
    if j == 0
        j = sitePerUnit;
    end
    switch dir
        case 'down'
            r = i;
            c = j;
        case 'right'
            r = j;
            c = i;
    end
end