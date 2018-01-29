function [volfrac,rho_2] = updatedensity_sensitivity_fsolve(Lambda,bn,P,vn,rho_min)
ne = size(bn,1);
rho_2 = zeros(ne,1);
for ie = 1:ne
    rho_target(ie) =(P*bn(ie)/(Lambda*vn(ie)))^(1/(P+1));
    rho_2(ie)= min([max([rho_target rho_min]) 1]);
end
volfrac = sum(rho_2)/ne;


