function transformout=splitmatch(pos1,pos2,dxv,maxd,transform,maxlocsused,show)
% function mytform=splitmatch(pos1,pos2,dx,dy,maxd,show,maxlocsused,transformationType)


if nargin<7||isempty(show)
    show=0;
end

if nargin<6||isempty(maxlocsused);
    maxlocsused=2000;
end

if nargin<5 || isempty(maxd)
    maxd=250;
end
if nargin<3 || isempty(dxv)
    dx=0;dy=0;
else
    dx=dxv(1);dy=dxv(2);
end


[iAa,iBa]=matchlocsall(pos1,pos2,dx,dy,maxd,maxlocsused);

disp([num2str(length(iAa)) ' of ' num2str(length(pos1.x))])

switch transform.type
    case {'lwm','polynomial'}
        mytform = fitgeotrans([double(pos1.x(iAa)) double(pos1.y(iAa))], [double(pos2.x(iBa)) double(pos2.y(iBa))],transform.type,transform.parameter);
    otherwise
        mytform = fitgeotrans([double(pos1.x(iAa)) double(pos1.y(iAa))], [double(pos2.x(iBa)) double(pos2.y(iBa))],transform.type);
end
% mytform = cp2tform([double(pos1.x(iAa)) double(pos1.y(iAa))], [double(pos2.x(iBa)) double(pos2.y(iBa))],'projective');


% mytform = cp2tform(double(pos1(iAa,2:3)), double(pos2(iBa,2:3)),'polynomial',3);
% mytform = cp2tform(double(pos1(iAa,2:3)), double(pos2(iBa,2:3)),'piecewise linear');
% mytform = cp2tform(double(pos1(iAa,2:3)), double(pos2(iBa,2:3)),'lwm');
transformout=transform;
transformout.T=mytform;

if show
%    [xa, ya]=tforminv(mytform,double(pos2(iBa,2)),double(pos2(iBa,3)));
%    tf.T=mytform;

   [xa, ya]=myapplytform((pos2.x(iBa)),(pos2.y(iBa)),transformout); 
   dx=xa-pos1.x(iAa);
   dy=ya-pos1.y(iAa);
%    figure(88)
   dscatter(dx,dy)
   title({['number of anchor points: ' num2str(length(iBa))],['dx= ' num2str(std(dx),3) ' nm, dy= ' num2str(std(dy),3) ' nm']});
 
end