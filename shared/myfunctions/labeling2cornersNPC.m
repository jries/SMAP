function labeling2cornersNPC(corners,rings,plot)
%calculate number of points in NPC (8x4) in dependence on labeling
%efficiency

plabel=0:0.1:1;

n=0:corners;


%probability that 1 corner is dark


% probability that n corners are dark
pnd=zeros(corners+1,length(plabel));
pnb=pnd;

for p=1:length(plabel)
    p1d=binopdf(0,rings,plabel(p));
    p1b=1-p1d;
    pnd(:,p)=binopdf(n,corners,p1d);
    %ncorers bright
    pnb(:,p)=binopdf(n,corners,p1b);
    
end

if plot
bar(plabel,pnb','stacked','DisplayName','pnbs')
title([num2str(corners) ' corners x ' num2str(rings) ' rings. number of bright corners'])
legend(num2str((0:corners)'),'Location','westoutside')
xlabel('labeling probability, single protein')
ylabel('probability')
xlim([-.05 1.05]);
ylim([-.05 1.05]);
end