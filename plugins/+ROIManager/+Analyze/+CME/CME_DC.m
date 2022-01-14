global se

figure(33);
hold off
for k=1%length(se.sites)
    eval=se.sites(k).evaluation.CME3DDSpherefit_DC;
    r1n=eval.pos1.r/eval.spherefit(1);
    r2n=eval.pos2.r/eval.spherefit(1);
    p1n=eval.pos1.p;
    p2n=eval.pos2.p;
    [x1,y1]=pol2cart(p1n,r1n);
    [x2,y2]=pol2cart(p2n,r2n);
    plot(x1,y1,'b.',x2,y2,'r.')
    hold on
end