function [z, bw] = getKernelMatrix(centralized, matrixSize, sVrot, varargin) 
    [meshx, meshy] = meshgrid(0:matrixSize(1), 0:matrixSize(2)); % make full grid
    if centralized
        shx = matrixSize(1)/2; shy = matrixSize(2)/2;
    else
        shx = 0; shy = 0;
    end
    
    q = [meshx(:)-shx, meshy(:)-shy]; % convert the grid to positions
    
    if length(varargin)==0
        [sVrotx,xy,bw]=ksdensity(sVrot, q); % get kernel density estimation
    else
        [sVrotx,xy,bw]=ksdensity(sVrot, q, varargin{1,1}, varargin{1,2}); % get kernel density estimation    
    end
    
    % converted into an image
    idx = sub2ind(matrixSize(end:-1:1)+1, xy(:,2)+1+shx, xy(:,1)+1+shy); 
    z = zeros(matrixSize(end:-1:1)+1);
    z(idx) = sVrotx;
end