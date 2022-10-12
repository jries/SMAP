function [z_phi,phicorr,phical] = z_from_phi_JR( z_ast, phi, k, z0 )
phical=(z_ast-z0)*2*k;
dphi=-(phi-phical)/pi/2;
period=round((dphi));
z_phi=period*pi/k+phi/2/k;
phicorr=mod(phi+z0*2*k,2*pi);
% figure(88);plot(z_ast,phi,'.',z_ast,phical,z_ast,period,'+',z_ast,dphi,'*')
end

