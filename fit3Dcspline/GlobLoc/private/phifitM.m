function [z_ang,z_phi_a,z_phi,phi] = phifitM(data4d,xf,yf,phi_s,phi_p,zT,M)
[rms, rmp] = iPALMast_findmom_givenM(data4d,xf,yf,M);
vec_amp = sqrt(rms.^2+rmp.^2);
rm1 = rms./vec_amp;
rm2 = rmp./vec_amp;
rm_complex = rm1-1i*rm2;
phi = angle(rm_complex);
z_phi = unwrap(phi).*zT./2./pi./1e3;

c=(rm1-1i*rm2);                                                            
z_angd=angle(c);
ang2=[];
parfor ii=1:numel(rms)
    [ang2(ii,:)]=iPALM_est_angle(double(rms(ii)),double(rmp(ii)),phi_s,phi_p,wrapToPi(double(z_angd(ii)-phi_s)));
end

z_ang = wrapToPi(ang2(:,2));
z_phi_a = unwrap(z_ang).*zT./2./pi./1e3;
