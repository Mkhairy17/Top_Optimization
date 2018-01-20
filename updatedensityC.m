function [volfrac,rho_2] = updatedensityC(Lambda,rho_1,R,a,b,nx,ny,P,strain,rho_min,Cm,F)
ne = size(strain,1);
rho_2 = zeros(ne,1);
Strain_energy = Calc_Strain_Energy(rho_1,strain,R,P,nx,ny,a,b,Cm,F);
for ie = 1:ne
    rho_target = (P*Strain_energy(ie)/Lambda)^(1/(P+1));
    rho_2(ie)= min([max([rho_target rho_min]) 1]);
end
volfrac = sum(rho_2)/ne;



