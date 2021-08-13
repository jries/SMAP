function out=simSplinePSF_call(Npixels,coeff,I,bg,cor)
if length(I)~=size(cor,1)
    I=zeros(size(cor,1),1,'single')+I(1);
end
if length(bg)~=size(cor,1)
    bg=zeros(size(cor,1),1,'single')+bg(1);
end
% corh=cor(:,[2 1 3]);
% corh(:,3)=-corh(:,3);
out=simSplinePSF_c(int16(Npixels),single(coeff),single(I),single(bg),single(cor));
% out=permute(out,[2 1 3]);