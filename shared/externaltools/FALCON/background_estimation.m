%% background estimation using iterative wavelet transform
function est_bg = background_estimation(imgs,th,dlevel,wavename,iter)
est_bg = zeros(size(imgs),'single');
imgs = max(imgs,0);
for N = 1: size(imgs,3)
    X = imgs(:,:,N);
    X_filt = X;
    
    for ii = 1:iter
        % wavelet transform
        [c,s] = wavedec2(X_filt,dlevel,wavename);
        cc = zeros(size(c));
        cc(1:s(1)*s(1)*1) = c(1:s(1)*s(1)*1);
        % inverse wavelet transform by only using low-freq components  
        X_new =  waverec2(cc,s,wavename);
        
        if th > 0
            % cut off values over current estimated background level.
            eps = sqrt(abs(X_filt))/2;
            ind = X>(X_new+eps);
            X_filt(ind) = X_new(ind)+eps(ind);

            % re-estimate background 
            [c,s] = wavedec2(X_filt,dlevel,wavename);
            cc = zeros(size(c));
            cc(1:s(1)*s(1)*1) = c(1:s(1)*s(1)*1);
            X_new =  waverec2(cc,s,wavename);
        end
    end
    est_bg(:,:,N) = X_new;
end
end
