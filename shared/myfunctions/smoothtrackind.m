function indstep=smoothtrackind(x,pin)
fac=2;
windowsize=ceil(length(x)/(max(x)-min(x))*pin.splitmergestep*fac);
time=(1:length(x))';

xf=runningWindowAnalysis(time,x,time,windowsize,pin.stepfunction);  
ind=1;
xh=xf(1);
indstep=0*xf;
for k=1:length(xf)
    if abs(xf(k)-xh)>pin.splitmergestep
        indstep(ind)=k;
        xh=xf(k);
        ind=ind+1;
    end

end
indstep(ind:end)=[];
end