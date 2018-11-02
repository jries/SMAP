function [xcenter,ycenter] = kernel_CenterofMass2D(sz,data)
tempx = single(0);
tempy = single(0);
tmpsum = single(0);
for ii = 0:sz-1
    for jj = 0:sz-1
        tempx = tempx+data(sz*jj+ii+1)*ii;
        tempy = tempy+data(sz*jj+ii+1)*jj;
        tmpsum = tmpsum +data(sz*jj+ii+1);
    end
end
xcenter = single(tempx/tmpsum);
ycenter = single(tempy/tmpsum);