function y=Hypoexp4fit(nt,k1,k2)
%fp: t1, t2,... tn, A 
    eps=0.000;
    ks=[k1 k2 k1+eps k2+eps];
%     A=1e6*0.2;
    y=Hypoexponentialpdf(nt,ks);
end