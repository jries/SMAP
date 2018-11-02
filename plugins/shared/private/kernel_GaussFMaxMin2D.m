function [MaxN,MinBG]=kernel_GaussFMaxMin2D(sz,sigma,data)
MaxN = single(0);
MinBG = single(10e10);
filteredpixel = single(0);
sum = single(0);
norm = single(1/2/sigma/sigma);

for kk = 0:sz-1
    for ll = 0:sz-1
        filteredpixel = single(0);
        sum = single(0);
        for ii = 0:sz-1
            for jj = 0:sz-1
                filteredpixel = filteredpixel+exp(-(ii-kk)^2*norm)*exp(-(ll-jj)^2*norm)*data(ii*sz+jj+1);
                sum = sum+exp(-(ii-kk)^2*norm)*exp(-(ll-jj)^2*norm);
            end
        end
        filteredpixel = filteredpixel/sum;
        MaxN = max(MaxN,filteredpixel);
        MinBG = min(MinBG, filteredpixel);
    end
end