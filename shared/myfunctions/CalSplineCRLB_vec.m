function [CRLB] =  CalSplineCRLB_vec(coeff, sz, Theta)
%input:
%coeff: spline coeffcients for PSF model
%sz: fit windowsize
%Theta: Nfit*5 array, collumn value [xpos,ypos,photon,bg,zpos]
%output:
%CRLB for the input Theta

spline_xsize = size(coeff, 1);
spline_ysize = size(coeff, 2);
spline_zsize = size(coeff, 3);

NV = 5; %number of parameters to optimize

off = floor(((spline_xsize+1)-sz)/2);

Nfits = size(Theta,1);
CRLB=zeros(Nfits,NV);
for tx = 1:Nfits
    xc=-1 * (Theta(tx,1) - sz / 2 + 0.5);
    yc = -1 * (Theta(tx,2) - sz / 2 + 0.5);
    zc = Theta(tx,5)-floor(Theta(tx,5));
    
    xstart=floor(xc);
    ystart=floor(yc);
    zstart=floor(Theta(tx,5));
    
    xc = xc-floor(xc);
    yc = yc-floor(yc);
    newTheta = Theta(tx,:);
%     [delta_f,delta_dxf,delta_ddxf1,delta_dyf,delta_ddyf1,delta_dzf,delta_ddzf1] = computeDelta3Dj_v2(single(xc),single(yc),single(zc));
% 
%     for ii = 1:sz
%         for jj = 1:sz
% 
%             
%             
%             [newDudt, model] =  kernel_DerivativeSpline_SMAP(ii+xstart+off,jj+ystart+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f,delta_dxf,delta_dyf,delta_dzf,coeff,newTheta);
%             
%             for l = 0:NV-1
%                 for m =l:NV-1
%                     
%                     M(l*NV+m+1) = M(l*NV+m+1)+newDudt(l+1)*newDudt(m+1)/model;
%                     M(m*NV+l+1) = M(l*NV+m+1);
%                 end
%             end
%             
%         end
%     end
    
    %again vectorized
    iiall=1:sz;
    jjall=1:sz;
    [delta_f,delta_dxf,delta_dyf,delta_dzf] = computeDelta3Dj_vec(single(xc),single(yc),single(zc));
    [newDudt, model] =  kernel_DerivativeSpline_SMAP_vec(iiall+xstart+off,jjall+ystart+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f,delta_dxf,delta_dyf,delta_dzf,coeff,newTheta);
    
    Mv = zeros(NV, NV);
    for l = 0:NV-1
        for m =l:NV-1

            Mv(l*NV+m+1) = sum(sum(newDudt(:,:,l+1).*newDudt(:,:,m+1)./model,1),2);%   Mv(l*NV+m+1)+newDudt(l+1)*newDudt(m+1)/model;
            Mv(m*NV+l+1) = Mv(l*NV+m+1);
        end
    end

    
    
    Minv = inv(Mv);
    for kk = 1:NV
        CRLB(tx,kk)=Minv(kk,kk);
    end
%     CRLBx(tx,:)=diag(Minv); %slower

end

% CRLB = reshape(CRLB,Nfits,NV);



