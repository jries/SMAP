function p=clusterfromlabeling(n,corners,rings,plabel)
%    n=1:corners;
% corners
% rings 
% plabel
   p1d=binopdf(0,rings,plabel);
   p1b=1-p1d;
   p=binopdf(n,corners,p1b);