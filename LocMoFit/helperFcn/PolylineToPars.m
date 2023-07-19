function out=PolylineToPars(input,angle)
%% Creates a structure with the initial parameters for continuoisLinearModel_PL_xyz
%input is a strusture created by PolylineToInit.m function
    %input must have fields called xnmR, ynmR, znmR of the same length
%out is a list of parameters that have the same number of control points as
%polyline

out=[];
ctrlPoints=size(input.xnmR,1);
coord=['x' 'y' 'z'];
if nargin<2
    angle=0;
end
for cp=1:ctrlPoints
    for i=coord
        field=[i 'nmR'];
        out(end+1).name=['m1.mPar.c' i num2str(cp)];
        out(end).value=input.(field)(cp);
    end
    if angle==1
        out(end+1).name=['m1.mPar.s' num2str(cp)];
        out(end).value=0;
    end
end
end
