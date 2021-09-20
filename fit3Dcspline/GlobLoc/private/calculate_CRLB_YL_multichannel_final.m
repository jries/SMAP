function [CRLB] =  calculate_CRLB_YL_multichannel_final(Nfits, coeff, sz, Theta,shared,dT)
% pi = single(3.141592);
spline_xsize = size(coeff, 1);
spline_ysize = size(coeff, 2);
spline_zsize = size(coeff, 3);

PSFSigma = single(1.5);
xc = single(0);
yc = single(0);
zc = single(0);
noChannels = 2;
sumShared = sum(shared);
NV = 5*noChannels-sumShared*(noChannels-1);

% NV = 6; %number of parameters to optimize

dudt = single(zeros(sz,sz,NV));


for i = 1:noChannels

    coeffall{i} = coeff(:,:,:,:,i);
end

off = floor(((spline_xsize+1)-sz)/2);
%
% Theta = zeros(NV,noChannels);

for tx = 1:Nfits
%     tx
    dS = dT(:,:,tx);
    M = zeros(NV, NV);
    Minv = zeros(NV, NV);
    for i = 1:noChannels
        
        %         Theta(1:2,i) = coordinates(1:2)
        xc(i)=single(-1 * (Theta(tx,1,i) - sz / 2 + 0.5));
        yc(i) = single(-1 * (Theta(tx,2,i) - sz / 2 + 0.5));
        zc(i) = single(Theta(tx,3,i)-floor(Theta(tx,3,i)));
        
        xstart(i)=floor(xc(i));
        ystart(i)=floor(yc(i));
        zstart(i)=floor(Theta(tx,3,i));
        
        xc(i) = xc(i)-floor(xc(i));
        yc(i) = yc(i)-floor(yc(i));
        
    end
    
    for i = 1: noChannels
        [delta_f(:,i),delta_dxf(:,i),delta_ddxf1,delta_dyf(:,i),delta_ddyf1,delta_dzf(:,i),delta_ddzf1] = computeDelta3Dj_v2(single(xc(i)),single(yc(i)),single(zc(i)));
        
        %      [delta_f2,delta_dxf2,delta_ddxf2,delta_dyf2,delta_ddyf2,delta_dzf2,delta_ddzf2] = computeDelta3Dj_v2(single(xc2),single(yc2),single(zc));
    end
    
    
    %     for i = 1: noChannels
    %         [delta_f(:,i),delta_dxf(:,i),delta_ddxf1,delta_dyf(:,i),delta_ddyf1,delta_dzf(:,i),delta_ddzf1] = computeDelta3Dj_v2(single(xc(i)),single(yc(i)),single(zc(i)));
    %
    %         %      [delta_f2,delta_dxf2,delta_ddxf2,delta_dyf2,delta_ddyf2,delta_dzf2,delta_ddzf2] = computeDelta3Dj_v2(single(xc2),single(yc2),single(zc));
    %     end
    %
    %
    %
    %     xc = single(-1*(Theta(tx, 1) - sz/2+0.5));
    %     yc = single(-1*(Theta(tx, 2) - sz/2+0.5));
    %     zc = single(Theta(tx, 5)-floor(Theta(tx, 5)));
    %     zstart = single(floor(Theta(tx, 5)));
    %
    %     xstart = floor(xc);
    %     xc = xc-floor(xc);
    %
    %     ystart = floor(yc);
    %     yc = yc-floor(yc);
    %     off = floor(((spline_xsize+1)-sz)/2);
    %     [delta_f,delta_dxf,delta_ddxf,delta_dyf,delta_ddyf,delta_dzf,delta_ddzf]=computeDelta3Dj_v2(single(xc),single(yc),single(zc));
    newDudt = zeros(5,noChannels);
    newDudtAll = zeros(5*noChannels*noChannels,1);
    for ii = 1:sz
        for jj = 1:sz
            
            for i = 1:noChannels
                
                %                 data2 = single(d_data2(sz*sz*(tx-1)+sz*jj+ii+1));
                
                delta_f1 = delta_f(:,i);
                delta_dxf1 = delta_dxf(:,i);
                delta_dyf1 = delta_dyf(:,i);
                delta_dzf1 = delta_dzf(:,i);
                %                 coeff1 = coeff(:,:,:,:,i);
                newTheta1 = Theta(tx,:,i);
                
                [newDudt1, model1] =  kernel_DerivativeSpline_v2_finalized(ii+xstart(i)+off,jj+ystart(i)+off,zstart(i),spline_xsize,spline_ysize,spline_zsize,delta_f1,delta_dxf1,delta_dyf1,delta_dzf1,coeffall{i},newTheta1,NV);
                newDudt(:,i)=newDudt1;
                model(i)=model1;
                %[newDudt(:,i), model(i)] =  kernel_DerivativeSpline_v2_finalized(ii+xstart(i)+off,jj+ystart(i)+off,zstart(i),spline_xsize,spline_ysize,spline_zsize,delta_f(:,i),delta_dxf(:,i),delta_dyf(:,i),delta_dzf(:,i),coeff(:,:,:,:,i),newTheta(:,i),NV);
                
                %                 [newDudt2, model2] =  kernel_DerivativeSpline_v2(ii+xstart2+off,jj+ystart2+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f2,delta_dxf2,delta_dyf2,delta_dzf2,coeff2,newTheta2,2);
            end
            %             [dudt, model] =  kernel_DerivativeSpline_v2(ii+xstart+off, jj+ystart+off, zstart, spline_xsize, spline_ysize, spline_zsize, delta_f, delta_dxf, delta_dyf, delta_dzf, PSF, squeeze(Theta(tx, :)), NV, phi0);
            
            n=1;
            for i = 1:5
                if shared(i) ==1
                    for j = 1:noChannels
                        newDudtAll(n+(j-1)*NV)=newDudt(i+(j-1)*5)*dS(noChannels*5+i+(j-1)*5);
                    end
                else
                    for j=1:noChannels
                        newDudtAll(n+j-1+(j-1)*NV)=newDudt(i+(j-1)*5);
                    end
                    n = n+j-1;
                end
                n=n+1;
                
            end
            
            
            
            for l = 0:NV-1
                for m =l:NV-1
                    for j = 1:noChannels
%                         M(l*NV+m+1) = M(l*NV+m+1)+newDudtAll(l+1,j)*newDudtAll(m+1,j)/model(j);
                        M(l*NV+m+1) = M(l*NV+m+1)+newDudtAll(l+1+(j-1)*NV)*newDudtAll(m+1+(j-1)*NV)/model(j);
                    end
                    M(m*NV+l+1) = M(l*NV+m+1);
                end
            end
            
        end
    end
    Minv = inv(M);
    for kk = 1:NV
        CRLB(Nfits*(kk-1)+tx)=Minv(kk,kk);
    end
    %             for kk = 1:1:length(phi0)
    %                 M(:, :, tx) = M(:, :, tx) + dudt(:, kk) * dudt(:, kk)' / model(kk);
    %             end
    %         end
    %     end
    %     Minv(:, :, tx) = inv(M(:, :, tx));
    %     for kk = length(squeeze(Theta(tx, :))):-1:1
    %         CRLB(tx,kk) = Minv(kk,kk,tx);
    %     end
end

CRLB = reshape(CRLB,Nfits,NV);

