function [Lambda] = initial_Lambda(rho_1,strain,P,a,b,Cm)
ne = size(strain,1);
for ie = 1:ne
   Strain_energy(ie) = rho_1(ie)^(2*P)*(strain(ie,:)*Cm*strain(ie,:)'*a*b); 
   Lambda0(ie) = (P/(rho_1(ie)^(P+1)))* Strain_energy(ie);
end
Lambda=mean(Lambda0);