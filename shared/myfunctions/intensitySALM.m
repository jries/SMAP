function [r,IS,IU]=intensitySALM(z,p)
if nargin==1
    n1=1.33;
    n2=1.78;
    lambda=680; %nm
    NA=1.70;
    NA0=1.33;

else
    n1=p.n1;
    n2=p.n2;
    lambda=p.lambda;
    NA=p.NA;
    NA0=p.NAmask;

end
lb=lambda/2/pi;
numax=sqrt(NA^2-n1^2)/lb;
wmax=n1/lb;

numin=0; 
wmin=0;
IUoff=0;
ISoff=0;
    
if NA0>n1
    numin=sqrt(NA0^2-n1^2)/lb;
    wmin=0;
    IUoff=integral(@dis,0,numin,'ArrayValued',true); %some part of supercritical emission ends up in undercritical channel
elseif NA0<n1
    numin=0;
    wmin=sqrt(n1^2-NA0^2)/lb;
    ISoff=integral(@diu,0,wmin); %some part of undercritical emission ends up in supercritical channel.
end

IS=integral(@dis,numin,numax,'ArrayValued',true)+ISoff;
IU=integral(@diu,wmin,wmax)+IUoff;
r=IS./IU;
% if limit
%     r(r>2)=2;
% end

    function is=dis(nu) %integration of supercritical part
        in = 2*(n1^2+n2^2)*nu*lb.*sqrt(n2^2-n1^2-nu.^2*lb^2).*(n1^1+nu.^2*lb^2);
        id=3*(n2^2-n1^2)*(n1^4+(n1^2+n2^2)*nu.^2*lb^2);
        is=in./id.*exp(-2*nu*z);
    end

    function iu=diu(w) %integration of undercritical part
        Q=sqrt(n2^2-n1^2+w.^2*lb^2);
        iu=2*Q*lb.*w/3.*(1./(w*lb+Q).^2+n1^2*n2^2./(n1^2*Q+n2^2*w*lb).^2);
    end
end