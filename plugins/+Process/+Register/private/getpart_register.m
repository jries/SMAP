function [ipart,m]=getpart_register(part,x,y,sim,midp)

switch part
    case 1 %all
        ipart=true(length(x),1);
        m.x=[0 sim(1)];
        m.y=[0 sim(2)];
    case 2% left
        ipart=x<midp;
        m.x=[0 midp];
        m.y=[0 sim(2)];
    case 3 %right
        ipart=x>=midp;
        m.x=[midp sim(1)];
        m.y=[0 sim(2)];
    case 4% up 
        ipart=y<midp;
        m.x=[0 sim(1)];
        m.y=[0 midp];
    case 5%down
        ipart=y>=midp;
        m.x=[0 sim(1)];
        m.y=[midp sim(2)];
end