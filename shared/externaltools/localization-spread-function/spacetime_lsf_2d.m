function [avglsf, lsfs] = spacetime_lsf_2d(x,y,t,spacewin,timewin,tau_edges,NmVal)
arguments
    x           (1,:)   double
    y           (1,:)   double
    t           (1,:)   double
    spacewin    (1,1)   struct {spacewin_isvalid}
    timewin     (:,2)   double {timewin_isvalid}
    tau_edges         (1,:)   double
    NmVal.XYpos    (1,:)   double = -50:4:50
    NmVal.How       (1,:)   string = 'actual'
end
[gs,~,~,Ns,~] = spacetime_acor_xy(x,y,t,spacewin,timewin,NmVal.XYpos,'TauEdges',tau_edges, 'How', NmVal.How);

lsfs = gs(:,:,1:end-1) - gs(:,:,end);
lsfs = lsfs ./ sum(sum(lsfs));

weights = sum(sum(Ns,1),2); % all pairs in given tau bin

avglsf = sum(cat(3,lsfs,lsfs(:,:,end)).*weights,3)/sum(weights,3);
end
