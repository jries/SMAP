function [P,CRLB,LogL]=mleFit_LM_HD_SE(imagestack,iterations,coeff,otherlocsROI,intipara)
% [P CRLB LL]=CPUmleFit_LM(single(output),5,iteration,single(coeff),isEMCCD,silent,otherlocsROI, initpara);
[P,CRLB,LogL]=CPUmleFit_HD(single(imagestack),5, iterations,single(coeff),0,1,single(otherlocsROI),single(intipara));