function ig=withinrec(x,y,crec,rimpix)
crec=single(crec);
if nargin==3
    rimpix=0;
end
rimpix=single(rimpix);
%crec 


ig=x>crec(3)-rimpix&x<crec(4)+rimpix&y>crec(1)-rimpix&y<crec(2)+rimpix;
% ig=ig&x<crec(4)+rimpix;
% ig=ig&y>crec(1)-rimpix;
% ig=ig&y<crec(2)+rimpix;
