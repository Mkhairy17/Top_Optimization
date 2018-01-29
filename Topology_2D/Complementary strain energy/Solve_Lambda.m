function [rho_new,Lambda_new] = Solve_Lambda(Lambda_old,Volume_Fraction_constraint,rho_old,P,rho_min,nx,ny,vn,bn)
Volume_Fraction_old = updatedensity_sensitivity(Lambda_old,bn,P,vn,rho_min);
while 1
drho_dLambda = -rho_old / ((P+1)*Lambda_old);
deta_dLambda = sum(drho_dLambda)/(nx*ny);
Lambda_new = Lambda_old - ( Volume_Fraction_old - Volume_Fraction_constraint)/(deta_dLambda);
[Volume_Fraction_new,rho_new] = updatedensity_sensitivity(Lambda_new,bn,P,vn,rho_min);
if norm(Volume_Fraction_new - Volume_Fraction_constraint,'inf') <= 1e-3
break;
end
rho_old=rho_new;
Lambda_old=Lambda_new;
Volume_Fraction_old = Volume_Fraction_new;
end