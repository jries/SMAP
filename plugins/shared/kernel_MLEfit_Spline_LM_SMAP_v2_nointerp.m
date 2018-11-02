function [P,CRLB,LogL] =  kernel_MLEfit_Spline_LM_SMAP_v2_nointerp(d_data,coeff,sz,iterations)
% update computeDelta3Dj and fAt3Dj. Add kernel_derivativeSpline
% make it capable of handling coeff with random size using Spline3D_v2
%work with odd number of ROI for calibration
CRLB=[];LogL=[];
pi = single(3.141592);
spline_xsize = size(coeff,1);
spline_ysize = size(coeff,2);
spline_zsize = size(coeff,3);

PSFSigma = single(1.5);
xc = single(0);
yc = single(0);
zc = single(0);
NV = 5;

dudt = single(zeros(sz,sz,NV));

newTheta = single(zeros(NV,1));
oldTheta = newTheta;

M = zeros(NV*NV,1);
Minv = zeros(NV*NV,1);

Nfits = size(d_data,3);

RUNNING = 0;
CONVERGED = 1;
CHOLERRER = 2;
BADPEAK = 3;
tolerance=1e-6;

for tx = 1:Nfits
    
    %newpeak
    [newTheta(1),newTheta(2)]=kernel_CenterofMass2D(sz,single(d_data(sz*sz*(tx-1)+1:sz*sz*(tx))));
    [Nmax,newTheta(4)] = kernel_GaussFMaxMin2D(sz,PSFSigma,single(d_data(sz*sz*(tx-1)+1:sz*sz*(tx))));
    newTheta(3)=max(0,(Nmax-newTheta(4))*2*pi*PSFSigma*PSFSigma);
    newTheta(5) = spline_zsize/2;
    
    newLambda = single(1);
    newStatus = RUNNING;
    newSign = zeros(NV,1);
    newUpdate = zeros(NV,1);
    newClamp = [1 1 max(100, newTheta(3)*50) 10000 2];
    %resetFit
    newFit = single(zeros(sz,sz));
    fit = single(zeros(sz,sz));
    
    %updateFitValues3D
%     xc = single(-2*(newTheta(1) - sz/2+0.5));
%     yc = single(-2*(newTheta(2) - sz/2+0.5));
    xc = single(-1*(newTheta(1) - sz/2+0.5));
    yc = single(-1*(newTheta(2) - sz/2+0.5));
    
    zc = single(newTheta(5)-floor(newTheta(5)));
    off = floor(((spline_xsize+1)-sz)/2);
    
    xstart = floor(xc);
    xc = xc-floor(xc);
    
    ystart = floor(yc);
    yc = yc-floor(yc);
    
    
    newErr = 0;
    jacobian = single(zeros(NV,1));
    hessian = single(zeros(NV,NV));
    zstart = single(floor(newTheta(5)));
    [delta_f,delta_dxf,delta_ddxf,delta_dyf,delta_ddyf,delta_dzf,delta_ddzf]=computeDelta3Dj_v2(single(xc),single(yc),single(zc));
    for ii = 0:sz-1
        for jj = 0:sz-1
            
            data = single(d_data(sz*sz*(tx-1)+sz*jj+ii+1));
            
            
            [newDudt,model] =  kernel_DerivativeSpline_v2(ii+xstart+off,jj+ystart+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f,delta_dxf,delta_dyf,delta_dzf,coeff,newTheta);
            
            
            
            if data>0
                newErr = newErr +2*((model-data)-data*log(model/data));
            elseif data ==0
                newErr = newErr + 2*model;
            end
            
            
            t1 = single(1-data/model);
            for l = 1:NV
                jacobian(l) = jacobian(l)+t1*newDudt(l);
            end
            t2 = data/model^2;
            
            for l = 0:NV-1
                for m =l:NV-1
                    hessian(l*NV+m+1) = hessian(l*NV+m+1)+t2*newDudt(l+1)*newDudt(m+1);
                    hessian(m*NV+l+1) = hessian(l*NV+m+1);
                end
            end
            
        end
    end
    %         dipshow(newDudt1);
    %     newErr = err;
    oldErr = 1e13;
    %addPeak
    
    %copyFitData
    
    oldClamp = newClamp;
    oldLambda = newLambda;
    oldSign = newSign;
    oldTheta = newTheta;
    
    for kk =1:iterations
        
        
        if abs((newErr-oldErr)/newErr)<tolerance
            %                 newStatus = CONVERGED;
            break;
        else
            if newErr>1.5*oldErr
                %subtractPeak
                fit = fit-newFit;
                %copyfitdata
                
                
                newClamp = oldClamp;
                newLambda = oldLambda;
                
                newSign = oldSign;
                
                newTheta = oldTheta;
                newErr = oldErr;
                
                
                newLambda = 10*newLambda;
                %addpeak
                
                
            elseif newErr<oldErr
                if newLambda >1
                    newLambda = newLambda*0.8;
                elseif newLambda <1
                    newLambda = 1;
                end
            end
            
            
            for i = 0:NV-1
                hessian(i*NV+i+1) = hessian(i*NV+i+1)*newLambda;
            end
            
            
            
            
            [L U info] = kernel_cholesky(hessian,NV);
            if info ==0
                newUpdate = kernel_luEvaluate(L,U,jacobian,NV);
                fit = fit -newFit;
                
                %copyFitData
                
                oldClamp = newClamp;
                oldLambda = newLambda;
                
                oldSign = newSign;
                
                oldTheta = newTheta;
                oldErr = newErr;
                
                
                %updatePeakParameters
                for ll =1:NV
                    if newSign(ll)~=0
                        if newSign(ll)==1&&newUpdate(ll)<0
                            newClamp(ll)=newClamp(ll)*0.5;
                        elseif newSign(ll)==-1&&newUpdate(ll)>0
                            newClamp(ll)=newClamp(ll)*0.5;
                        end
                        
                    end
                    
                    if newUpdate(ll)>0
                        newSign(ll)=1;
                    else
                        newSign(ll)=-1;
                    end
                    newTheta(ll) = newTheta(ll)-newUpdate(ll)/(1+abs(newUpdate(ll)/newClamp(ll)));
%                     update(kk,ll) = newUpdate(ll)/(1+abs(newUpdate(ll)/newClamp(ll)));
                end
%                 update(kk,NV+1)= oldErr;
                
                %                 newTheta(1) = max(newTheta(1),(sz-1)/2-2);
                %                 newTheta(1) = min(newTheta(1),(sz-1)/2+2);
                %                 newTheta(2) = max(newTheta(2),(sz-1)/2-2);
                %                 newTheta(2) = min(newTheta(2),(sz-1)/2+2);
                newTheta(3) = max(newTheta(3),1);
                %                 newTheta(4) = max(newTheta(4),0.01);
                newTheta(5) = max(newTheta(5),0);
                newTheta(5) = min(newTheta(5),spline_zsize);
                
                %updateFitValues3D
                %                 xc = single(-2*((newTheta(1) - sz/2)+0.5));
                %                 yc = single(-2*((newTheta(2) - sz/2)+0.5));
                %                 zc = single(newTheta(5)-floor(newTheta(5)));
                
                xc = single(-1*(newTheta(1) - sz/2+0.5));
                yc = single(-1*(newTheta(2) - sz/2+0.5));
                zc = single(newTheta(5)-floor(newTheta(5)));
                
                xstart = floor(xc);
                xc = xc-floor(xc);
                
                ystart = floor(yc);
                yc = yc-floor(yc);
                
                newErr = 0;
                jacobian = single(zeros(NV,1));
                hessian = single(zeros(NV,NV));
                zstart = single(floor(newTheta(5)));
                [delta_f,delta_dxf,delta_ddxf,delta_dyf,delta_ddyf,delta_dzf,delta_ddzf]=computeDelta3Dj_v2(single(xc),single(yc),single(zc));
                for ii = 0:sz-1
                    for jj = 0:sz-1
                        [newDudt,model] =  kernel_DerivativeSpline_v2(ii+xstart+off,jj+ystart+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f,delta_dxf,delta_dyf,delta_dzf,coeff,newTheta);
                        data = single(d_data(sz*sz*(tx-1)+sz*jj+ii+1));
                        
                        if data>0
                            newErr = newErr +2*((model-data)-data*log(model/data));
                        elseif data ==0
                            newErr = newErr + 2*model;
                        end
                        
                        t1 = single(1-data/model);
                        for l = 1:NV
                            jacobian(l) = jacobian(l)+t1*newDudt(l);
                        end
                        t2 = data/model^2;
                        
                        for l = 0:NV-1
                            for m =l:NV-1
                                hessian(l*NV+m+1) = hessian(l*NV+m+1)+t2*newDudt(l+1)*newDudt(m+1);
                                hessian(m*NV+l+1) = hessian(l*NV+m+1);
                            end
                        end
                        
                    end
                end
                %                                 dipshow(newDudt1);
            else
                newLambda = 10*newLambda;
                disp('CHOLERRER')
            end
            %                 kk
            
            
        end
        
    end
    iteration = kk;
    for kk = 1:NV
        P(Nfits*(kk-1)+tx) = newTheta(kk);
    end
    P(Nfits*(NV)+tx) = iteration;
end
P = reshape(P,Nfits,NV+1);








































