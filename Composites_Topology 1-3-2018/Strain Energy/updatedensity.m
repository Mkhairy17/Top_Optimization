function [volfrac,rho_new] = updatedensity(Lambda,rho_old,strain,P,rho_min,Cm,a,B) 
zeta = 0.05;
eta = 1;
ne = size(strain,1);
rho_new = zeros(ne,1);
for ie = 1:ne
    Strain_energy=a*B*strain(ie,:)*Cm*strain(ie,:)';
    b =  P*rho_old(ne)^(P-1)*Strain_energy/Lambda;
    rho_target = b^eta*rho_old(ie);
    rho_u = min([(1+zeta)*rho_old(ie) 1]);
    rho_l = max([(1-zeta)*rho_old(ie) rho_min])  ;  
    rho_new(ie) = max(min(rho_target,rho_u),rho_l);
end

volfrac = sum(rho_new)/ne;

% end

