function dat=fitimagesimple(imstack,p)
f=figure;
wffile='settings/workflows/onlyfit.mat';
wf=interfaces.Workflow(f);
wf.makeGui;
wf.load(wffile);
loader=wf.module('Arrayfeeder');

fitter=wf.module('MLE_GPU_Yiming');
collector=wf.module('LocCollector');
peakfinder=wf.module('PeakFinder');

loader.imstack=single(imstack);
% pfit=struct('roisperfit',100,'iterations',30);
if isfield(p,'mindistance')
    p.use_mindistance=true;
end
peakfinder.setGuiParameters(p);
fitter.setGuiParameters(p);
wf.run;

dat=collector.locdata.loc;
close(f);
end