global se
sites=se.sites;
stepWidthAll = [];
stallTimeAll = [];
growthTimeAll = [];
avgRateAll = [];
idxAll = [];
k = 1;
dxn=diff(sites(k).evaluation.BALM_fibril_growth.xn);
dxn=dxn(1);
dfr=diff(sites(k).evaluation.BALM_fibril_growth.fr);
dfr=dfr(1);
dt=dfr*150/1000/60;
for k=1:length(sites)
    name=sites(k).name;
    if ~isequal(regexp(name, '\+$'),[])
    try
        stepWidth=sites(k).evaluation.boundaryFinder.stepWidth*dxn;
        stallTime=sites(k).evaluation.boundaryFinder.stallTime*dt;
        growthTime=sites(k).evaluation.boundaryFinder.growthTime*dt;
        avgRate=sites(k).evaluation.boundaryFinder.avgRate*dxn/dt;
        
    catch err
        continue
    end
    stepWidthAll = [stepWidthAll; stepWidth(2:(end-1))];
    stallTimeAll = [stallTimeAll; stallTime(2:(end-1))];
    growthTimeAll = [growthTimeAll; growthTime(2:(end-1))];
    avgRateAll = [avgRateAll; avgRate];
    idxAll = [idxAll, k];
    
    end
end
growthRateAll = stepWidthAll./growthTimeAll;

figure(500)
subplot(2,3,1);
histogram(stallTimeAll, 'Binwidth', dt*2);
title('Stall time')
xlabel('min')

subplot(2,3,2);
histogram(growthTimeAll, 'Binwidth', dt);
title('Growth time')
xlabel('min')

subplot(2,3,3);
histogram(stepWidthAll, 'Binwidth', dxn);
title('Step width')
xlabel('nm')

subplot(2,3,4);
histogram(growthRateAll, 'Binwidth', dxn/dt/2);
title('Growth rate')
xlabel('nm/min')

subplot(2,3,5);
histogram(avgRateAll, 'Binwidth', dxn/dt/20);
title('Average rate')
xlabel('nm/min')

subplot(2,3,6);
scatter(growthTimeAll, stepWidthAll, 'jitter','on', 'jitterAmount', 1);
title('Growth time vs step width')
xlabel('Growth time (min)')
ylabel('step width (nm)')

corr(growthTimeAll, stepWidthAll)