function likelihood = gaussDist(refLocs,y,x,z,ysigmao,xsigmao,zsigmao, varargin)
    p = inputParser;
    addParameter(p, 'sigFactor',1);
    addParameter(p, 'rotMat',[]);
    addParameter(p, 'rotMode','old');
    addParameter(p, 'distMode','ordinary');
    addParameter(p, 'distCutoff',3.5);
    addParameter(p, 'sumIntensity',[]);
    
    parse(p, varargin{:})
    % For one dimension: exp(-((x-mu/sigma)^2)/2)/(2*pi)^(1/2)*sigma
    f = p.Results.sigFactor;          % scale factor of the 3 sigma
    rotMat = p.Results.rotMat;        % rotation matrix
    distMode = p.Results.distMode; 
    rotMode = p.Results.rotMode;
    distCutoff = p.Results.distCutoff;
    sumIntensity = p.Results.sumIntensity;
    if isempty(sumIntensity)
        sumIntensity = 1-(1-mvncdf(repelem(((p.distCutoff^2)/3)^(1/2),1,3),[0 0 0],diag([1 1 1])))*2;
    end
    
    switch rotMode
        case 'old'
            xsigma = xsigmao.*f;ysigma = ysigmao.*f;zsigma = zsigmao.*f;

            % the template bellow is for multiple sets of parameters at the same
            % time
            templateAll = zeros([1 size(refLocs.x,2) size(refLocs.x,1)]);   % zero filled vector
            template.x = templateAll;
            template.y = templateAll;
            template.z = templateAll;
            template.x(1,:,:)=refLocs.x';
            template.y(1,:,:)=refLocs.y';
            template.z(1,:,:)=refLocs.z';
            refLocs = template;
            
            switch distMode
                % Fast mode saves time by only evaluating nearby locs
                % around the model points.
                case 'fast'
                    dx = (x - refLocs.x)./xsigma;
                    dy = (y - refLocs.y)./ysigma;
                    dz = (z - refLocs.z)./zsigma;
%                     lSkip = dx>2*xsigma|dy>2*ysigma|dz>2*zsigma;
%                     dx(lSkip) = inf;
                    termR2 = dx.^2 + dy.^2 + dz.^2;
%                     termR2 = reshape(termR2, size(termR2,[1 3]));
                    svol = ((2*pi)^(-3/2))./(xsigma.*ysigma.*zsigma);
                    
                    expDist = zeros(size(termR2));
                    l = termR2<distCutoff^2;
%                     expDist(l) = exp(-termR2(l)/2)./sumIntensity;
                    expDist(l) = exp(-termR2(l)/2);
                    
%                     dim = size(expDist);
%                     expDist = reshape(expDist, [dim(1) 1 dim(2)]);
                case 'ordinary'
                    dx = x./xsigma - refLocs.x./xsigma;
                    dy = y./ysigma - refLocs.y./ysigma;
                    dz = z./zsigma - refLocs.z./zsigma;
                    termR2 = dx.^2 + dy.^2 + dz.^2;
                    svol = ((2*pi)^(-3/2))./(xsigma.*ysigma.*zsigma);
                    expDist = exp(-termR2/2);

        %             (x - refLocs.x).^2+(y - refLocs.y).^2+(z - refLocs.z).^2
        %             likelihood = likelihood';
            end
            likelihood = svol.*sum(expDist,3);
        case 'new'
            %% rotation of SIGMA (covariance matrix)
            % rotMat
            rotMat = rotMat(1:3,1:3,:);
            rotMat_originalSize = size(rotMat,3);
            
            if 0
                % this part is only for validation
                rotMat = diag([1,1,1]);
                rotMat = repmat(rotMat, [1 1 rotMat_originalSize]);
            end
            
            % SIGMA
            SIGMA = zeros([3 3 length(ysigmao)]);
            SIGMA_originalSize = size(SIGMA,3);
            refCopy_originalSize = size(refLocs.y,1);
            refDiff_originalSize = size(refLocs.y,2);
            SIGMA(1,1,:)=xsigmao*f.^2; SIGMA(2,2,:)=ysigmao*f.^2; SIGMA(3,3,:)=zsigmao*f.^2;
            SIGMA = repmat(SIGMA,[1 1 max([rotMat_originalSize refDiff_originalSize])]);
            
            % SIGMArot
            if  rotMat_originalSize>1
                rotMat = repelem(rotMat,1, 1, SIGMA_originalSize);
            else
                rotMat = repelem(rotMat,1, 1, SIGMA_originalSize*refDiff_originalSize);
            end
            rotMat_inv = permute(rotMat,[2 1 3]);
            currentTerm = ndfun('mult',SIGMA, rotMat_inv);
            SIGMArot = ndfun('mult',rotMat, currentTerm);
            ref = [refLocs.x refLocs.y refLocs.z];
            locs = [x y z];
            % 
            dist = locs'-permute(ref, [2 3 1]);
            dist = permute(dist, [1 4 2 3]);
%             dist(1,:) = x(:); dist(2,:) = y(:); dist(3,:) = z(:);
            % add info of ref
            SIGMArot_inv = multinv(SIGMArot);
            SIGMA_det = SIGMArot(1,1,:).*SIGMArot(2,2,:).*SIGMArot(3,3,:)+...
             SIGMArot(1,2,:).*SIGMArot(2,3,:).*SIGMArot(3,1,:)+...
             SIGMArot(1,3,:).*SIGMArot(2,1,:).*SIGMArot(3,2,:)-...
             SIGMArot(1,3,:).*SIGMArot(2,2,:).*SIGMArot(3,1,:)-...
             SIGMArot(1,1,:).*SIGMArot(2,3,:).*SIGMArot(3,2,:)-...
             SIGMArot(1,2,:).*SIGMArot(2,1,:).*SIGMArot(3,3,:);
            SIGMA_det_sqrt = SIGMA_det.^-0.5;
            SIGMA_det_sqrt = squeeze(SIGMA_det_sqrt);
%             SIGMArot_inv = permute(SIGMArot_inv, [1 2 4 3]);
%             t1 = SIGMArot_inv.*dist;
            distT = permute(dist,[2 1 3 4]);
%             currentTerm = ndfun('mult',SIGMArot_inv, dist);
            t2 = repmat(SIGMArot_inv);
            t1 = SIGMArot_inv.*dist;
            t1 = sum(t1,1);
            t1 = sum(distT.*t1,2);
            t1 = squeeze(t1);
            t1 = exp(-t1./2);
            currentTerm = (2*pi).^-(3/2)*SIGMA_det_sqrt.*t1;
            likelihood = sum(currentTerm,2);
    end
end
