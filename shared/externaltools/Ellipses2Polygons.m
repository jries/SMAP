function Polygons=Ellipses2Polygons(Ellipses,L)

%Input
%Ellipses: (NEllipses*5) array. each column represent: 
           %(:,1): radius of ellipses at direction 1
           %(:,2): radius of ellipses at direction 2
           %(:,3): x-cooridante of centroid
           %(:,4): y-coordinate of centroid
           %(:,5): Inclination angle
%L: Approximate segement length of the polygons

%Output
%Polygons: a cell array where Polygons{i} is (n*2) polygons vertices array of the ith ellipse 
           
m=size(Ellipses,1);

for i=1:1:m    
    
r1=Ellipses(i,1);
r2=Ellipses(i,2);
Cx=Ellipses(i,3);
Cy=Ellipses(i,4); 
alpha=Ellipses(i,5);

P=2*pi*sqrt((r1^2+r2^2)/2);
N=P/L;
N=4*round(N/4); 
if N==0; N=4; end

th=zeros(N/4,1);
for j=1:1:N/4
   A=(r1+r2)*pi/4;
   rts=roots([(r1-r2)/pi r2 -4*A/N*j]);
   th(j)=rts(2);
end
the=[th ; pi-th(end-1:-1:1) ; pi ; th+pi ; 2*pi-th(end-1:-1:1) ; 2*pi];

x1=r1*cos(the); y1=r2*sin(the);

x2=x1*cos(alpha)-y1*sin(alpha);
y2=x1*sin(alpha)+y1*cos(alpha);

x2=x2+Cx;
y2=y2+Cy;

v=[x2 y2];

Polygons{i}=v;

end

end