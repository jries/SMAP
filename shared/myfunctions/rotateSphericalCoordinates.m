function [elr,azimutr]=rotateSphericalCoordinates(eli,azimut,axis,alpha)
switch axis
    case 1 %x
        R=[cos(alpha/2) -1i*sin(alpha/2)
            -1i*sin(alpha/2) cos(alpha/2)];
    case 2 %y
         R=[cos(alpha/2) -sin(alpha/2)
            sin(alpha/2) cos(alpha/2)];
        
    case 3 %z
%          R=[exp(-1i*alpha/2) 0
%             0 exp(1i*alpha/2)];
        elr=eli;
        azimutr=azimut+alpha;
        return
end
elir=reshape(eli,1,numel(eli));
azimutr=reshape(azimut,1,numel(azimut));

[elrm,azimutrm]=rotatei(elir,azimutr);

elr=reshape(elrm,size(eli));
azimutr=reshape(azimutrm,size(azimut));
% [elr,azimutr]=arrayfun(@rotatei,eli,azimut);

% sum(elrd-elr);
% sum(azimutrd-azimutr)

function [elr,azimutr]=rotatei(eli,azimut)
theta2=pi/4-eli/2;
z1=[cos(theta2); exp(1i*azimut).*sin(theta2)];




z2=R*z1;
elr=2*atan2(abs(z2(2,:)),abs(z2(1,:)));
azimutr=angle(z2(2,:))-angle(z2(1,:));

elr=pi/2-elr;
end

end