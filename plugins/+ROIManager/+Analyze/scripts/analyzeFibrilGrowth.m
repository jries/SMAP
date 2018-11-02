global se
frametime=se.locData.files.file(1).info.timediff/1000; %in s
sites=se.sites;
dfall=[];
vall=[];
lenall=[];
timeall=[];
vaverage=[];
for k=1:length(sites)
    try
        polyh=sites(k).evaluation.BALM_fibril_growth.poly;
        dpol=polyh(2:end,:)-polyh(1:end-1,:);
    catch err
        continue
    end
    l=size(dpol,1);
    dfall(end+1:end+l)=dpol(:,2)*frametime;
     lenall(end+1:end+l)=dpol(:,1)/1000;
    vall(end+1:end+l)=dpol(:,1)./dpol(:,2)/frametime/1000; %in um/s
    timeall(end+1:end+l)=polyh(1:end-1,2)*frametime;
    
    %average growth speed
    try
        dpol=polyh(end-1,:)-polyh(2,:);
        vaverage(k)=dpol(:,1)./dpol(:,2)/frametime/1000;
    catch err
        continue
    end
end

maxwaitelongation=0.03; %in um
stop=lenall<maxwaitelongation;


figure(99)
subplot(2,3,1);
histogram(dfall(stop),0:1000:max(dfall(stop)));
title('wait time')
xlabel('wait time s')
subplot(2,3,2);
histogram(dfall(~stop),0:1000:max(dfall(~stop)));
title('growth time')
xlabel('growth time s')


subplot(2,3,3);
histogram(lenall(~stop),0:.010:max(lenall(~stop)));
title('elongationstep')
xlabel('elongationstep um')

subplot(2,3,4);
histogram(vall(~stop),0:1e-5:min(max(vall(~stop)),.01));
title('elongation speed')
xlabel('elongation speed um/s')

subplot(2,3,5);
histogram(vaverage*1000*60,0:1e-1:min(quantile(vaverage(~isinf(vaverage)),.99),.01)*1000*60);
title('average elongation speed')
xlabel('average speed nm/min')

figure(100);
plot(timeall(stop),dfall(stop),'+')



