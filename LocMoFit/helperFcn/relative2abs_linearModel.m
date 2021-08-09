function parOut = relative2abs_linearModel(varargin)
    nArgin = length(varargin);
    if rem(nArgin,3)==0
        x0 = varargin{1};
        y0 = varargin{2};
        z0 = varargin{3};
        
        nRestSet = nArgin/3-1;
        for k = nRestSet:-1:1
            dist(k) = varargin{3*(k)+1}; % ->
            rotAzi(k) = varargin{3*(k)+2}; % ->
            rotEle(k) = varargin{3*(k)+3}; % ->
        end
        rotAzi = cumsum(rotAzi',1);
        rotEle = cumsum(rotEle',1);
        rotAzi = rotAzi.*pi./180;
        rotEle = rotEle.*pi./180;

        [diffX,diffY,diffZ] = sph2cart(rotAzi,rotEle,dist');
        %% 191128
        diffX = [x0;diffX];
        diffY = [y0;diffY];
        diffZ = [z0;diffZ];
        ctrlX = cumsum(diffX,1);
        ctrlY = cumsum(diffY,1);
        ctrlZ = cumsum(diffZ,1);
        
        for k = (nRestSet:-1:1)+1
            parOut.(['x' num2str(k-1)]) = ctrlX(k);
            parOut.(['y' num2str(k-1)]) = ctrlY(k);
            parOut.(['z' num2str(k-1)]) = ctrlZ(k);
        end
    else
        warning('Input control points might not be complete.')
    end
end