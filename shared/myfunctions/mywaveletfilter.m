function out=mywaveletfilter(in,level,refine,cutoff)
if nargin<4
    cutoff=false;
end
if nargin<3
    refine=false;
end
persistent wf1 wf2
if isempty(wf1)
    wf1=load('near_sym_a.mat');
    wf2=load('qshift_a.mat');
end
if nargin<2
    level=3;
end
if nargin<3
    refine=false;
end
% ino=in;
% if cutoff
%     co=myquantilefast(in(:),.99);
%     in(in>co)=co;
% end
if refine
    level2=level*2;
    [L,S]=dtwavexfm2_L(in,level2,wf1,wf2);
    bg=dtwaveifm2_L(L,level2,wf1,wf2,S);
    if numel(bg)~=numel(in)
        s=size(in);
        bg=bg(1:s(1),1:s(2));
    end
    in=in-bg;co=3*std(in(:)); in(in>co)=co;   
end
[L,S]=dtwavexfm2_L(in,level,wf1,wf2);
out=dtwaveifm2_L(L,level,wf1,wf2,S);
if numel(out)~=numel(in)
    s=size(in);
    out=out(1:s(1),1:s(2));
end
if refine
    out=out+bg;
end


if 0
H0=3/8;
H1=1/4;
H2=1/16;
g1=[H2,H1,H0,H1,H2];
g2=[H2,0,H1,0,H0,0,H1,0,H2];

V1=conv2(conv2(in,g1','same'),g1,'same');
V2=conv2(conv2(V1,g2','same'),g2,'same');
out=im-V2;
end