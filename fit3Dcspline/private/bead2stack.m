
function [stack,filenumber,dT]=bead2stack(beads)
ss=size(beads(1).stack.image);
if length(ss)==3
    stack=zeros(ss(1),ss(2),ss(3),length(beads));
    for k=length(beads):-1:1
        stack(:,:,:,k)=beads(k).stack.image;
        filenumber(k)=beads.filenumber;
    end
elseif length(ss)==4
    stack=zeros(ss(1),ss(2),ss(3),length(beads),ss(4));
    numpar=6;
    dT=zeros(numpar,ss(4),ss(3),length(beads));
    for k=length(beads):-1:1
        stack(:,:,:,k,:)=beads(k).stack.image;
        filenumber(k)=beads.filenumber;
        for zz=1:ss(3)
            dT(1:2,:,zz,k)=squeeze(beads(k).shiftxy(1,[2 1],:));
        end
    end    
end
end