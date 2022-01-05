
function [P,CRLB,LogL] = psfloc_4Pi(obj,rois,dx,dy,shared)
data = rois;
data = data+obj.BGoffset;
% IratioL = zeros(1,1,1,4);
% IratioL(1,1,1,:) = obj.Iratio;
% bxsz = obj.Boxsize;
numlocs = size(rois,3);
% data = data./repmat(IratioL,[bxsz,bxsz,numlocs,1]);

zstart = ones(numlocs,1,'single')*[obj.Initz];
p0 = ones(numlocs,1,'single')*single(obj.InitPhase);
sharedA = repmat(shared',[1 numlocs]);
phi0 = obj.Phi0;
phi0A = repmat(phi0',[1 numlocs]);

LogL = ones(numlocs,1,'single')-1e10;
sz0 = size(zstart);
sp0 = size(p0);
iterations = obj.iterations;

P = zeros(numlocs,7+3*(6-sum(shared)));
CRLB = zeros(numlocs,6+3*(6-sum(shared)));

dTAll = cat(1,reshape(dx',1,4,[]),reshape(dy',1,4,[]),zeros(1,4,numlocs),zeros(1,4,numlocs),repmat(obj.Dz,1,1,numlocs),repmat(obj.Dphi,1,1,numlocs));


for p=1:sp0(2)
    for k=1:sz0(2)
        [Ph,CRLBh, LogLh] = CPUmleFit_LM_4Pi_nIAB(single(data),uint32(sharedA),iterations,single(obj.IABall),single(dTAll),single(phi0A),single(zstart(:,k)'),single(p0(:,p)'));
        indbetter=LogLh-LogL>1e-4; %copy only everything if LogLh increases by more than rounding error.
        P(indbetter,:)=Ph(indbetter,:);
        CRLB(indbetter,:)=CRLBh(indbetter,:);
        LogL(indbetter)=LogLh(indbetter);
    end
end
    

clear CPUmleFit_LM_4Pi_nIAB 
  
    

