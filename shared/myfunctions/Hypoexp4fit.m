function y=Hypoexp4fit(nt,k1,k2,A)
%fp: t1, t2,... tn, A 
    eps=0.01;
    ks=[k1 k2 k1+eps k2+eps];
%     A=1e6*0.2;
A=1;
    y=A*Hypoexponentialpdf(nt,ks);
end