function outim= invertwhite(in,p)
if nargin<2
    p=0;
end
if isinteger(in)
    maxn=2^(8*prod(st.size)/st.bytes);
    in=double(in)/maxn;
end
eps=1e-12;
% inc=in/max(in(:));
inc=in;
inc(inc>1)=1;
inten=max(eps,max(inc,[],3));
% inten2=min(max(eps,sum(inc,3)),2);
% inten=(inten+inten2)/2;
color=in./inten;
colori=max(1-color,0);
graygood=(sum(colori,3)>0.5);
colori=colori.*graygood+~graygood;
% colori=colori./sum(colori,3);
colori(isnan(colori))=1;
outim=1-inten.*colori;
% outim=(1+(inten.*in-inten)*0.9)-inten.*colori;
outim=outim.*(1-inten*p);

if isinteger(in)
    outim=outim*maxn;
end
outim=cast(outim,'like',in);