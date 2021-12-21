%analysis of MINFLUX tracks
sites=g.locData.SE.sites;

stepsize=[]; %define the variables you want to read out
steptime=[];

pluginname='StepsMINFLUX';
fieldname1='steps';


for k=1:length(sites)
    if ~isfield(sites(k).evaluation,pluginname) || ~isfield(sites(k).evaluation.(pluginname),fieldname1) % if no steps are found, look at next site
        continue
    end
    sh=sites(k).evaluation.(pluginname).(fieldname1);
    stepsize(end+1:end+length(sh.stepsize))=sh.stepsize; %add values from current site to the list
    steptime(end+1:end+length(sh.dwelltime))=sh.dwelltime;
end

figure(88) %do the plotting
subplot(2,2,1)
ds=1;
n=round(min(stepsize)):ds:max(stepsize);
histogram(stepsize,n)
xlabel('stepsize(nm)')

dt=5;
n=0:dt:max(steptime);
subplot(2,2,2)
histogram(steptime,n)
xlabel('step time (ms)')