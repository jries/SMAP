function [gpus,gpurec]=parsegpustathttp(gpustat)
fn=fieldnames(gpustat);
for k=1:length(fn)
    gpus.mem(k)=gpustat.(fn{k}).memory_total-gpustat.(fn{k}).memory_util;
    gpus.load(k)=gpustat.(fn{k}).load;
    gpus.name{k}=strrep(fn{k},'_',':');
end

[mmax, ind]=max(gpus.mem);
gpurec=gpus.name{ind};
end