function [xpos,fp]=mypeakfit(x,y)
fp=polyfit(x,y,2);
xpos=-fp(2)/2/fp(1);
end