function indstep=findstepsMINFLUX(x,pin)
%p. splitmerge (T), splitmergestep([]), stepfunction (mean),
p.fitmean=contains(pin.stepfunction,'mean');
p.splitmerge=pin.splitmerge;
p.splitmergestep=pin.splitmergestep;
p.estimatenoise=false;
p.bootstraprepeats=0;
try
steps=AutoStepfinderRies(x,p);
catch err
    disp('AutoStepFinder did not work')
    steps=[];
end
if isempty(steps) %try splitmerge
%     disp('no steps found')
    steps.properties.StepSize=10;
    steps.properties.IndexStep=ceil(length(x)/2);
    steps.properties.LevelBefore=x(1);
    steps.properties.LevelAfter=x(end);

end
indstep=[1 ;steps.properties.IndexStep ];      

if p.splitmerge
    stepsize=p.splitmergestep;
    if isempty(stepsize)
        stepsize=median(steps.properties.StepSize);
    end
    [indstep,stepvalue]=splitmergefit(x,stepsize,p,steps);
end
end


function [istepfit,svalfit]=splitmergefit(x,stepsize,p,steps)
if p.fitmean
    mfun=@mysimplemean;
else
    mfun=@median;
end

istep=[1 ; steps.properties.IndexStep];
sval=[steps.properties.LevelBefore; steps.properties.LevelAfter(end)];
    svalfit=sval;
    istepfit=istep;
for s=1:10
    [istepfit, svalfit]=mergesplit(istepfit,svalfit,stepsize);
    % recalculate sval based on x and istep2
    [svalfit, istepfit]=fitstepind(x,istepfit,mfun);
   
end
end







function [istep, sval]=mergesplit(istep,sval,stepsize)
stepv=diff(sval);
step2=stepv(1:end-1)+stepv(2:end);
sstep=find(abs(step2)<stepsize*1.4 & abs(stepv(1:end-1))<stepsize*.7);
k=1;
sstepc=sstep;
while k<length(sstepc)
    ind=find(sstepc==sstepc(k)+1);
    if ~isempty(ind)
        sstepc(ind)=[];
        %k=k+1;
    end
    k=k+1;
end
istep(sstepc+2)=round((istep(sstepc+1)+istep(sstepc+2))/2);
istep(sstepc+1)=[];
sval(sstepc+1)=[];
bigstep=find(abs(diff(sval))>stepsize*1.4)+1; %later pass on min / max step size

for k=1:length(bigstep)
    [sval,istep]=insertstep(sval,istep,bigstep(k));
    bigstep=bigstep+1; % as inserted, correct indices
end
end

function [sval,istep]=insertstep(sval,istep,insertind)
%     if insertind>=length(istep)
%         return
%     end
    istep=[istep(1:insertind) ;istep(insertind) ;istep(insertind+1:end)];
    sval=[sval(1:insertind-1) ;(sval(insertind-1)+sval(insertind))/2 ;sval(insertind:end)];  
    istep(insertind)=istep(insertind)-1;
    istep(insertind+1)=istep(insertind+1)+1;
end

% function [sval,istep]=removestep(sval,istep,insertind)
% if insertind+1<=length(istep)
%     istep(insertind+1)=round((istep(insertind)+istep(insertind+1))/2);
% end
%     istep(insertind)=[];
%     sval(insertind)=[];
% end
