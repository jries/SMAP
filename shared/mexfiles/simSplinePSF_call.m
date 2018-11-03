function out=simSplinePSF_call(Npixels,coeff,I,bg,cor)
out=simSplinePSF_c(int16(Npixels),single(coeff),single(I),single(bg),single(cor));