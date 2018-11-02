function ind=mywithin(x,cx,y,cy)
%c: left, width
ind=(x>cx(1)&x<cx(1)+cx(2));

if nargin==4
    indg2=(y>cy(1)&y<cy(1)+cy(2));
    ind=ind&indg2;
end
