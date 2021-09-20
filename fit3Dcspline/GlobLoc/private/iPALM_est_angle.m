% for v14 iPALMast_analysisv14scmos.m

function [x]=iPALM_est_angle(rms,rmp,phi1,phi2,ang0)
a0=0.5;
options = optimset('Display','off');
x=fsolve(@fun,[a0 ang0],options);

    function f=fun(x)
        f(1)=abs(x(1))*cos(x(2)+phi1)-rms;
        f(2)=abs(x(1))*cos(x(2)+phi2)-rmp;
    end
end
% r=rms./rmp;
% tan_theta=(r.*cos(phi2)+cos(phi1))./(r.*sin(phi2)+sin(phi1));
% ang=atan(tan_theta);