function lut=mymakelut(lutname)
if nargin==0
    lut={'red hot','green cold','cyan cold','jet', 'red','green','blue','magenta','cyan','yellow','orange','purple','bgy','bry','gray','gray inv','hsv','parula','parula black','kgy'};
else
switch lutname
    case 'red hot'
        lut=hot(256);
    case 'red'
        lut=zeros(256,3);
        lut(:,1)=(0:255)/255;
    case 'green'
        lut=zeros(256,3);
        lut(:,2)=(0:255)/255;
    case 'blue'
        lut=zeros(256,3);
        lut(:,3)=(0:255)/255;    
    case 'cyan'
        lut=zeros(256,3);
        lut(:,3)=(0:255)/255; 
        lut(:,2)=(0:255)/255; 
    case 'magenta'
        lut=zeros(256,3);
        lut(:,1)=(0:255)/255; 
        lut(:,3)=(0:255)/255; 
   case 'purple'
        lut=zeros(256,3);
        lut(:,1)=(0:255)/255; 
        lut(:,3)=(0:255)/412;
    case 'yellow'
        lut=zeros(256,3);
        lut(:,1)=(0:255)/255; 
        lut(:,2)=(0:255)/255; 
    case 'orange'
        lut=zeros(256,3);
        lut(:,1)=(0:255)/255; 
        lut(:,2)=(0:255)/384; 
    case 'jet'
        lut=jet(256);
    case 'green cold'
        lut=usercolormap([0 0 0],[0 1 0],[0 1 1],[1 1 1]);
    case 'cyan cold'
        lut=usercolormap([0 0 0],[0 .2 1],[0 0.5 1],[0 1 1],[.5 1 1],[1 1 1]);
    case 'kgy'
        lut=usercolormap([0 0 0],[.5 .25 0],[1 .5 0],[1 1 0],[1 1 1]);
    case 'bgy'
        lut=usercolormap([0 0 1],[0 1 0],[1 1 0]);
    case 'bry'
        lut=usercolormap([0 0 1],[1 0 0],[1 1 0]);
    case 'gray'
        lut=gray(256);
    case 'gray inv'
        lut=gray(256);
        lut=lut(end:-1:1,:);
    case 'hsv'
        lut=hsv(256);
    case 'parula'
        lut=parula(256);
    case 'parula black'
        c=parula(14);
        lut=usercolormap([0 0 0],c(1,:),c(2,:),c(3,:),c(4,:),c(5,:),c(6,:),c(7,:),c(8,:),c(9,:),c(10,:),c(11,:),c(12,:),c(13,:),c(14,:));
%         lut=usercolormap([0 0 0],[0 0.2 1],[0.0 .5 .7],[0.2 .9 .5],[.5 .9 0.2],[.7 .7 0.0],[1 .7 0.0],[1 1 0],[1 1 .2],[1 1 .7]);

end
end