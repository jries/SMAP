function plotsquares(x,y,hl,attribute,ploth)
for k=1:length(x)
     plotrect(ploth, [x(k)-hl,y(k)-hl,x(k)+hl,y(k)+hl],attribute); 
end
end